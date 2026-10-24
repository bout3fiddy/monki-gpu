// SPDX-License-Identifier: MIT
// Adapted from MONKI 1.0.1, Copyright (c) 2026 Victor J.H. Trees.
// Source map: mc.f90; mc_routines_general_mod.f90; mc_routines_cartesian_mod.f90.
// Forward I/Q/U transport, periodic Cartesian grid, flat Lambertian surface.
// The original Fortran prepares the optical tables.
// Platform types, math, address spaces and atomics are supplied by the adapter.
MONKI_CONSTANT float PI = 3.14159265358979323846f;
MONKI_CONSTANT float EPS = 1.0e-16f;
MONKI_CONSTANT uint MAX_LAYERS = 35;

struct Scene {
    int nx, ny, nz, nmat, angles, lut_size;
    bool paths;
    float xunit, yunit, depol, theta0, phi0, phi, aerosol_albedo;
    float fixed_views[MONKI_ANGLES];
    MONKI_GLOBAL const float *zunit, *views, *albedo, *beta, *cloud_fraction, *ssalb;
    MONKI_GLOBAL const float *cosines, *cdf, *f11, *f22, *f33, *f12;
    MONKI_GLOBAL const float *height_above, *tau_above;
};

MONKI_DEVICE Scene scene_tables(MONKI_GLOBAL const float *data) {
    Scene s;
    s.nx = MONKI_NX;
    s.ny = MONKI_NY;
    s.nz = MONKI_NZ;
    s.nmat = MONKI_NMAT;
    s.angles = MONKI_ANGLES;
    s.lut_size = MONKI_LUT_SIZE;
    s.paths = MONKI_PATHS != 0;
    s.xunit = data[7];
    s.yunit = data[8];
    s.depol = data[9];
    s.theta0 = MONKI_THETA0;
    s.phi0 = MONKI_PHI0;
    s.phi = MONKI_PHI;
    s.aerosol_albedo = data[13];

    s.zunit = data + 16;
    s.views = s.zunit + s.nz;
    s.fixed_views[0] = MONKI_VIEW_0;
#if MONKI_ANGLES > 1
    s.fixed_views[1] = MONKI_VIEW_1;
#endif
#if MONKI_ANGLES > 2
    s.fixed_views[2] = MONKI_VIEW_2;
#endif
    s.albedo = s.views + s.angles;
    s.beta = s.albedo + s.nx * s.ny;
    s.cloud_fraction = s.beta + s.nx * s.ny * s.nz;
    s.ssalb = s.cloud_fraction + s.nx * s.ny * s.nz;
    s.cosines = s.ssalb + s.nz;
    s.cdf = s.cosines + s.lut_size;
    s.f11 = s.cdf + s.lut_size;
    s.f22 = s.f11 + s.lut_size;
    s.f33 = s.f22 + s.lut_size;
    s.f12 = s.f33 + s.lut_size;
    s.height_above = s.f12 + s.lut_size;
    s.tau_above = s.height_above + s.nz;
    return s;
}

// Philox4x32-10: streams keyed by seed, round, photon identity, and draw block.
// Scheduling/chunk size never changes a photon's stream; this intentionally does
// not imitate gfortran RANDOM_NUMBER. See Random123's published Philox constants.
struct RandomStream {
    uint photon, round, seed, block, lane;
    uint4 words;
};

MONKI_DEVICE uint4 philox(uint4 counter, uint2 key) {
    for (uint i = 0; i < 10; ++i) {
        MONKI_UINT64 p0 = MONKI_UINT64(0xD2511F53u) * counter.x;
        MONKI_UINT64 p1 = MONKI_UINT64(0xCD9E8D57u) * counter.z;
        counter = uint4(uint(p1 >> 32) ^ counter.y ^ key.x, uint(p1),
                        uint(p0 >> 32) ^ counter.w ^ key.y, uint(p0));
        key += uint2(0x9E3779B9u, 0xBB67AE85u);
    }
    return counter;
}

MONKI_DEVICE float random_uniform(MONKI_THREAD RandomStream &rng) {
    if (rng.lane == 4) {
        rng.words = philox(uint4(rng.photon, rng.round, rng.block++, 0), uint2(rng.seed, 0));
        rng.lane = 0;
    }

    // Exactly representable half steps on an open (0,1) float32 grid: no log(0)
    // and no rounded-to-one endpoint in inverse-CDF sampling.
    return (float(rng.words[rng.lane++] >> 9) + 0.5f) * 0x1p-23f;
}

MONKI_DEVICE void add_tally(MONKI_GLOBAL atomic_float *target, float value) {
    atomic_fetch_add_explicit(target, value, memory_order_relaxed);
}

MONKI_DEVICE int bracket(MONKI_GLOBAL const float *axis, int count, float value) {
    int low = 0, high = count - 1;
    while (high - low > 1) {
        int middle = (low + high) / 2;
        if (axis[middle] <= value)
            low = middle;
        else
            high = middle;
    }
    return low;
}

// ray_scat_matrix / mie_scat_matrix: return F11,F22,F33,F12, including /4
// normalization (MONKI radiance is normalized assuming incident E0 = pi).
MONKI_DEVICE float4 scattering_elements(MONKI_THREAD const Scene &s, bool gas, float cosine) {
    cosine = clamp(cosine, -1.0f, 1.0f);
    float4 f;
    if (gas) {
        float delta = (1.0f - s.depol) / (1.0f + s.depol / 2.0f);
        float diagonal = delta * 0.75f * (1.0f + cosine * cosine);
        f = float4(diagonal + 1.0f - delta, diagonal, delta * 1.5f * cosine,
                   -delta * 0.75f * (1.0f - cosine * cosine));
    } else {
        int low = bracket(s.cosines, s.lut_size, cosine);
        float fraction = (cosine - s.cosines[low]) / (s.cosines[low + 1] - s.cosines[low]);
        f = mix(float4(s.f11[low], s.f22[low], s.f33[low], s.f12[low]),
                float4(s.f11[low + 1], s.f22[low + 1], s.f33[low + 1], s.f12[low + 1]), fraction);
    }
    if (s.nmat == 1)
        f.w = 0;
    return f / 4.0f;
}

// stokes_rot_mat: retain Victor's coplanar, vertical, and general cases.
// The epsilon for angular degeneracy is float32-scale, unlike the FP64 source.
MONKI_DEVICE float4 rotation(float cosine, float uo, float u1, float phio, float phi1) {
    float dphi = abs(phio - phi1);
    float remainder = fmod(dphi, PI);
    if (min(remainder, PI - remainder) < 2e-7f || (abs(uo) == 1.0f && abs(u1) == 1.0f))
        return float4(1, 0, 1, 0);

    float sino = sqrt(max(0.0f, 1.0f - uo * uo));
    float sin1 = sqrt(max(0.0f, 1.0f - u1 * u1));
    float2 wophi(cos(phio), sin(phio)), w1phi(cos(phi1), sin(phi1));
    float2 wo = sino * wophi, w1 = sin1 * w1phi;
    float cosi1, cosi2, q;
    if (abs(u1) == 1.0f) {
        cosi1 = copysign(1.0f, u1);
        cosi2 = -cosi1 * dot(w1phi, wophi);
        q = dot(float2(wo.y, -wo.x), w1phi);
    } else if (abs(uo) == 1.0f) {
        cosi2 = copysign(1.0f, uo);
        cosi1 = -cosi2 * dot(w1phi, wophi);
        q = dot(float2(wophi.y, -wophi.x), w1);
    } else {
        cosi1 = (u1 - uo * cosine) / sqrt(max(EPS, (1 - uo * uo) * (1 - cosine * cosine)));
        cosi2 = (uo - u1 * cosine) / sqrt(max(EPS, (1 - u1 * u1) * (1 - cosine * cosine)));
        q = dot(float2(wo.y, -wo.x), w1);
    }

    cosi1 = clamp(cosi1, -1.0f, 1.0f);
    cosi2 = clamp(cosi2, -1.0f, 1.0f);
    float sini1 = copysign(sqrt(max(0.0f, 1 - cosi1 * cosi1)), q);
    float sini2 = copysign(sqrt(max(0.0f, 1 - cosi2 * cosi2)), q);
    return float4(cosi1 * cosi1 - sini1 * sini1, 2 * sini1 * cosi1, cosi2 * cosi2 - sini2 * sini2,
                  2 * sini2 * cosi2);
}

// matmul_stokes_rot followed by P*S, with the original I/Q/U sign convention.
MONKI_DEVICE float3 scatter_stokes(MONKI_THREAD const Scene &s, bool gas, float cosine, float uo,
                                   float u1, float phio, float phi1, float3 stokes) {
    float4 f = scattering_elements(s, gas, cosine);
    float4 r = rotation(cosine, uo, u1, phio, phi1);
    float rotated_q = r.x * stokes.y - r.y * stokes.z;
    float rotated_u = r.y * stokes.y + r.x * stokes.z;
    float intensity = f.x * stokes.x + f.w * rotated_q;
    float q = f.w * stokes.x + f.y * rotated_q;
    float u = f.z * rotated_u;
    return float3(intensity, r.z * q - r.w * u, r.w * q + r.z * u);
}

// rayl_scatter_pol / mie_scatter_pol_bound, forward sign = -1.
// A rejected sampling limit is a reported failure, never a biased fallback.
MONKI_DEVICE bool sample_direction(MONKI_THREAD const Scene &s, bool gas, float3 stokes, float uo,
                                   MONKI_THREAD RandomStream &rng, MONKI_THREAD float &u1,
                                   MONKI_THREAD float &dphi, MONKI_THREAD float &sampled_cosine) {
    for (uint attempt = 0; attempt < 4096; ++attempt) {
        float cosine;
        if (gas)
            cosine = 2 * random_uniform(rng) - 1;
        else {
            float draw = random_uniform(rng);
            int low = bracket(s.cdf, s.lut_size, draw);
            float fraction = (draw - s.cdf[low]) / (s.cdf[low + 1] - s.cdf[low]);
            cosine = mix(s.cosines[low], s.cosines[low + 1], fraction);
        }

        float sigma = 2 * PI * random_uniform(rng) - PI;
        float4 f = scattering_elements(s, gas, cosine);
        float pdf = f.x + f.w * (stokes.y * cos(2 * sigma) - stokes.z * sin(2 * sigma));
        float bound = gas ? scattering_elements(s, true, 1).x : f.x + abs(f.w);
        if (random_uniform(rng) * bound > pdf)
            continue;

        sampled_cosine = cosine;
        float sin_in = sqrt(max(EPS, 1 - uo * uo));
        float sin_sca = sqrt(max(EPS, 1 - cosine * cosine));
        u1 = clamp(uo * cosine + sin_in * sin_sca * cos(sigma), -1.0f, 1.0f);
        float fraction = (cosine - uo * u1) / (sin_in * sqrt(max(EPS, 1 - u1 * u1)));
        dphi = -acos(clamp(fraction, -1.0f, 1.0f)) * copysign(1.0f, sigma);
        return true;
    }
    return false;
}

struct Photon {
    int i, j, k;
    float3 remaining; // distances to the cell faces in the current direction
    float3 direction;
    float phi;
};

MONKI_DEVICE int cell_index(MONKI_THREAD const Scene &s, MONKI_THREAD const Photon &p) {
    return p.i + s.nx * (p.j + s.ny * p.k);
}

MONKI_DEVICE void change_direction(MONKI_THREAD const Scene &s, MONKI_THREAD Photon &p,
                                   float3 direction, float phi) {
    float3 size(s.xunit, s.yunit, s.zunit[p.k]);
    p.remaining =
        select(p.remaining, max(float3(0), size - p.remaining), p.direction * direction < 0);
    p.direction = direction;
    p.phi = phi;
}

MONKI_DEVICE float boundary_distance(MONKI_THREAD const Photon &p) {
#if MONKI_NX == 1 && MONKI_NY == 1
    // A periodic single-column atmosphere cannot change at an x/y face.
    // Skip those bookkeeping crossings while retaining the full slant distance.
    return p.remaining.z / max(abs(p.direction.z), EPS);
#else
    float3 distances = p.remaining / max(abs(p.direction), float3(EPS));
    return min(distances.x, min(distances.y, distances.z));
#endif
}

// ADVANCE_internal: preserve the x-before-y-before-z boundary tie convention.
MONKI_DEVICE float advance(MONKI_THREAD const Scene &s, MONKI_THREAD Photon &p) {
#if MONKI_NX == 1 && MONKI_NY == 1
    float distance = boundary_distance(p);
    p.k += p.direction.z >= 0 ? 1 : -1;
    if (p.k >= 0 && p.k < s.nz)
        p.remaining.z = s.zunit[p.k];
    return distance;
#else
    float3 distances = p.remaining / max(abs(p.direction), float3(EPS));
    float distance = min(distances.x, min(distances.y, distances.z));
    p.remaining = max(float3(0), p.remaining - abs(p.direction) * distance);
    if (distance == distances.x) {
        p.i = (p.i + (p.direction.x >= 0 ? 1 : -1) + s.nx) % s.nx;
        p.remaining.x = s.xunit;
    } else if (distance == distances.y) {
        p.j = (p.j + (p.direction.y >= 0 ? 1 : -1) + s.ny) % s.ny;
        p.remaining.y = s.yunit;
    } else {
        p.k += p.direction.z >= 0 ? 1 : -1;
        if (p.k >= 0 && p.k < s.nz)
            p.remaining.z = s.zunit[p.k];
    }
    return distance;
#endif
}

// Project a detector ray through the uniform upper atmosphere. This is the
// original Cartesian endpoint projection, with periodic wrapping before indexing.
MONKI_DEVICE int projected_column(int index, float remaining, float direction, float distance,
                                  float unit, int count) {
    float local = direction >= 0 ? unit - remaining : remaining;
    float position = fmod(float(index) * unit + local + direction * distance, unit * count);
    if (position < 0)
        position += unit * count;
    return min(int(position / unit), count - 1);
}

// calc_radiances_and_paths_pol: event estimation along the detector ray.
// Only the top-of-atmosphere output grid is supported; traversing to its boundary
// gives the same periodic pixel as the original analytical projection.
MONKI_DEVICE bool detector(MONKI_THREAD const Scene &s, MONKI_THREAD const Photon &event, int angle,
                           float3 probability, MONKI_THREAD const float *history,
                           MONKI_GLOBAL atomic_float *tallies, uint max_steps) {
    if (probability.x <= EPS)
        return true;
    Photon ray = event;
    float mu = cos(s.fixed_views[angle]);
    float3 view(sin(s.fixed_views[angle]) * cos(s.phi), sin(s.fixed_views[angle]) * sin(s.phi), mu);
    change_direction(s, ray, view, s.phi);
    float tau = 0;
    float lengths[MAX_LAYERS] = {};
    uint steps = 0;

    while (ray.k >= 0 && ray.k < s.nz - 1) {
        if (++steps > max_steps)
            return false;
        if (ray.k >= MONKI_UNIFORM_TOP) {
            // Every remaining layer has identical extinction across x/y. The
            // vertical optical-depth prefix therefore gives the exact slant
            // integral; retain each layer's length for the AMF estimator.
            float distance = (ray.remaining.z + s.height_above[ray.k]) / mu;
            tau += (ray.remaining.z * s.beta[cell_index(s, ray)] + s.tau_above[ray.k]) / mu;
            if (s.paths) {
                lengths[ray.k] += ray.remaining.z / mu;
                for (int k = ray.k + 1; k < s.nz - 1; ++k)
                    lengths[k] += s.zunit[k] / mu;
            }
            ray.i = projected_column(ray.i, ray.remaining.x, view.x, distance, s.xunit, s.nx);
            ray.j = projected_column(ray.j, ray.remaining.y, view.y, distance, s.yunit, s.ny);
            break;
        }
        int cell = cell_index(s, ray), layer = ray.k;
        float distance = advance(s, ray);
        tau += distance * s.beta[cell];
        if (s.paths)
            lengths[layer] += distance;
    }

    float3 contribution = probability * exp(-tau);
    int pixel = ray.i + s.nx * (ray.j + s.ny * angle);
    float geometry = abs(cos(s.theta0)) / mu;
    for (int component = 0; component < s.nmat; ++component)
        add_tally(tallies + s.nmat * pixel + component, contribution[component] * geometry);

    if (s.paths) {
        int stokes_size = s.nmat * s.nx * s.ny * s.angles;
        int path_size = s.nx * s.ny * s.nz * s.angles;
        for (int k = 0; k < s.nz; ++k) {
            int offset = ray.i + s.nx * (ray.j + s.ny * (k + s.nz * angle));
            add_tally(tallies + stokes_size + offset, contribution.x * (history[k] + lengths[k]));
        }
        add_tally(tallies + stokes_size + path_size + pixel, contribution.x);
    }
    return true;
}

MONKI_DEVICE void transport_photon(MONKI_GLOBAL const float *data,
                                   MONKI_GLOBAL atomic_float *tallies,
                                   MONKI_GLOBAL atomic_uint *failure,
                                   MONKI_CONSTANT_PTR uint *controls, uint lane) {
    if (lane >= controls[0])
        return;
    Scene s = scene_tables(data);
    // Stripe photons across independent output banks. Only the final reduction
    // changes; each photon keeps the same RNG stream and physical estimator.
    uint tally_count = s.nmat * s.nx * s.ny * s.angles;
    if (s.paths)
        tally_count += s.nx * s.ny * (s.nz + 1) * s.angles;
    uint identity = controls[1] + lane;
    uint round_offset = identity / controls[6];
    uint photon_id = identity % controls[6];
    // Grouped rounds share scheduling only: RNG keys and tally banks stay separate.
    tallies += (round_offset * uint(data[14]) + lane % uint(data[14])) * tally_count;
    RandomStream rng = {photon_id, controls[2] + round_offset, controls[3], 0, 4, uint4(0)};
    uint max_steps = controls[4];
    int column = photon_id % uint(s.nx * s.ny);
    Photon p;
    p.i = column % s.nx;
    p.j = column / s.nx;
    p.k = s.nz - 1;
    p.direction = float3(sin(s.theta0) * cos(s.phi0), sin(s.theta0) * sin(s.phi0), cos(s.theta0));
    p.phi = s.phi0;
    p.remaining =
        float3(s.xunit * random_uniform(rng), s.yunit * random_uniform(rng), s.zunit[p.k]);
    float3 stokes(1, 0, 0);
    float weight = 1;
    float history[MAX_LAYERS] = {};
    float remains = -log(random_uniform(rng));

    for (uint step = 0; step < max_steps; ++step) {
        float extinction = p.k < s.nz - 1 ? s.beta[cell_index(s, p)] : 0;
        float event_distance = extinction > EPS ? remains / extinction : INFINITY;
        bool surface = false;
        if (event_distance < boundary_distance(p)) {
            p.remaining = max(float3(0), p.remaining - abs(p.direction) * event_distance);
            if (s.paths)
                history[p.k] += event_distance;
        } else {
            int old_layer = p.k;
            float distance = advance(s, p);
            if (s.paths)
                history[old_layer] += distance;
            remains -= distance * extinction;
            if (p.k >= s.nz)
                return;
            if (p.k >= 0)
                continue;
            surface = true;
            p.k = 0;
            p.remaining.z = 0;
        }

        bool gas = true;
        if (surface)
            weight *= s.albedo[p.i + s.nx * p.j];
        else {
            gas = random_uniform(rng) >= s.cloud_fraction[cell_index(s, p)];
            weight *= gas ? s.ssalb[p.k] : s.aerosol_albedo;
        }
        if (weight < EPS)
            return;

        for (int angle = 0; angle < s.angles; ++angle) {
            float mu = cos(s.fixed_views[angle]);
            float3 probability;
            if (surface)
                probability = float3(stokes.x * mu * weight, 0, 0);
            else {
                float cosine =
                    clamp(p.direction.z * mu + sqrt(max(0.0f, 1 - p.direction.z * p.direction.z)) *
                                                   sin(s.fixed_views[angle]) * cos(p.phi - s.phi),
                          -1.0f, 1.0f);
                probability =
                    scatter_stokes(s, gas, cosine, p.direction.z, mu, p.phi, s.phi, stokes) *
                    weight;
            }
            if (!all(isfinite(probability))) {
                atomic_fetch_or_explicit(failure, 1u, memory_order_relaxed);
                return;
            }
            if (!detector(s, p, angle, probability, history, tallies, max_steps)) {
                atomic_fetch_or_explicit(failure, 2u, memory_order_relaxed);
                return;
            }
        }

        float mu, phi;
        if (surface) {
            mu = sqrt(random_uniform(rng));
            phi = 2 * PI * random_uniform(rng);
            stokes = float3(1, 0, 0);
        } else {
            float dphi, cosine;
            if (!sample_direction(s, gas, stokes, p.direction.z, rng, mu, dphi, cosine)) {
                atomic_fetch_or_explicit(failure, 4u, memory_order_relaxed);
                return;
            }
            phi = p.phi + dphi;
            // Use the sampled scattering cosine directly, as in mc.f90;
            // reconstructing it from float32 directions loses forward-peak precision.
            stokes = scatter_stokes(s, gas, cosine, p.direction.z, mu, p.phi, phi, stokes);
            stokes /= stokes.x;
            if (!all(isfinite(stokes)) || sqrt(stokes.y * stokes.y + stokes.z * stokes.z) > 1.01f) {
                atomic_fetch_or_explicit(failure, 8u, memory_order_relaxed);
                return;
            }
            remains = -log(random_uniform(rng));
        }

        float sine = sqrt(max(0.0f, 1 - mu * mu));
        change_direction(s, p, float3(cos(phi) * sine, sin(phi) * sine, mu), phi);
    }

    atomic_fetch_or_explicit(failure, 16u, memory_order_relaxed);
}
