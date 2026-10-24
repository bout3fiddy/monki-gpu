// CUDA adaptations for the small Metal vector/math vocabulary used by MONKI.
// Physics and RNG equations live in common/transport.h. All operations are FP32.
#include <cuda_runtime.h>
#define MONKI_DEVICE __device__ inline
#define MONKI_CONSTANT static constexpr
#define MONKI_GLOBAL
#define MONKI_THREAD
#define MONKI_CONSTANT_PTR const
#define INFINITY __int_as_float(0x7f800000)
using uint = unsigned int;
using MonkiUint64 = unsigned long long;
#define MONKI_UINT64 MonkiUint64
#define abs monki_abs
using atomic_float = float;
using atomic_uint = unsigned int;
constexpr int memory_order_relaxed = 0;

// Explicit vectors avoid relying on CUDA's plain floatN structs for operators.
struct Vec2 {
    float x, y;
    MONKI_DEVICE Vec2(float a, float b) : x(a), y(b) {}
};
struct Vec3 {
    float x, y, z;
    MONKI_DEVICE Vec3() = default;
    MONKI_DEVICE Vec3(float a) : x(a), y(a), z(a) {}
    MONKI_DEVICE Vec3(float a, float b, float c) : x(a), y(b), z(c) {}
    MONKI_DEVICE float &operator[](int i) {
        return i == 0 ? x : (i == 1 ? y : z);
    }
    MONKI_DEVICE float operator[](int i) const {
        return i == 0 ? x : (i == 1 ? y : z);
    }
    MONKI_DEVICE Vec3 &operator/=(float a) {
        x /= a;
        y /= a;
        z /= a;
        return *this;
    }
};
struct Vec4 {
    float x, y, z, w;
    MONKI_DEVICE Vec4() = default;
    MONKI_DEVICE Vec4(float a, float b, float c, float d) : x(a), y(b), z(c), w(d) {}
};
struct Words2 {
    uint x, y;
    MONKI_DEVICE Words2(uint a, uint b) : x(a), y(b) {}
    MONKI_DEVICE Words2 &operator+=(Words2 b) {
        x += b.x;
        y += b.y;
        return *this;
    }
};
struct Words4 {
    uint x, y, z, w;
    MONKI_DEVICE Words4(uint a) : x(a), y(a), z(a), w(a) {}
    MONKI_DEVICE Words4(uint a, uint b, uint c, uint d) : x(a), y(b), z(c), w(d) {}
    MONKI_DEVICE uint operator[](uint i) const {
        return i == 0 ? x : (i == 1 ? y : (i == 2 ? z : w));
    }
};
#define float2 Vec2
#define float3 Vec3
#define float4 Vec4
#define uint2 Words2
#define uint4 Words4
MONKI_DEVICE Vec2 operator*(float a, Vec2 b) {
    return Vec2(a * b.x, a * b.y);
}
MONKI_DEVICE Vec3 operator*(Vec3 a, Vec3 b) {
    return Vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}
MONKI_DEVICE Vec3 operator*(Vec3 a, float b) {
    return Vec3(a.x * b, a.y * b, a.z * b);
}
MONKI_DEVICE Vec3 operator-(Vec3 a, Vec3 b) {
    return Vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}
MONKI_DEVICE Vec3 operator/(Vec3 a, Vec3 b) {
    return Vec3(a.x / b.x, a.y / b.y, a.z / b.z);
}
MONKI_DEVICE Vec4 operator/(Vec4 a, float b) {
    return Vec4(a.x / b, a.y / b, a.z / b, a.w / b);
}
struct Mask3 {
    bool x, y, z;
};
MONKI_DEVICE Mask3 operator<(Vec3 a, float b) {
    return {a.x < b, a.y < b, a.z < b};
}
MONKI_DEVICE Vec3 select(Vec3 a, Vec3 b, Mask3 m) {
    return Vec3(m.x ? b.x : a.x, m.y ? b.y : a.y, m.z ? b.z : a.z);
}
MONKI_DEVICE float abs(float a) {
    return fabsf(a);
}
MONKI_DEVICE Vec3 abs(Vec3 a) {
    return Vec3(fabsf(a.x), fabsf(a.y), fabsf(a.z));
}
MONKI_DEVICE Vec3 max(Vec3 a, Vec3 b) {
    return Vec3(fmaxf(a.x, b.x), fmaxf(a.y, b.y), fmaxf(a.z, b.z));
}
MONKI_DEVICE float clamp(float a, float low, float high) {
    return fminf(fmaxf(a, low), high);
}
MONKI_DEVICE float dot(Vec2 a, Vec2 b) {
    return a.x * b.x + a.y * b.y;
}
MONKI_DEVICE float mix(float a, float b, float t) {
    return a + (b - a) * t;
}
MONKI_DEVICE Vec4 mix(Vec4 a, Vec4 b, float t) {
    return Vec4(mix(a.x, b.x, t), mix(a.y, b.y, t), mix(a.z, b.z, t), mix(a.w, b.w, t));
}
MONKI_DEVICE bool all(bool a) {
    return a;
}
MONKI_DEVICE bool isfinite(Vec3 a) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
}
MONKI_DEVICE void atomic_fetch_add_explicit(float *p, float v, int) {
    atomicAdd(p, v);
}
MONKI_DEVICE void atomic_fetch_or_explicit(uint *p, uint v, int) {
    atomicOr(p, v);
}
