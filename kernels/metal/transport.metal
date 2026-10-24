// Apple entry point; shared equations retain their original MONKI routine mapping.
#include <metal_stdlib>
using namespace metal;
#define MONKI_CONSTANT constant
#define MONKI_DEVICE inline
#define MONKI_UINT64 ulong
#define MONKI_GLOBAL device
#define MONKI_THREAD thread
#define MONKI_CONSTANT_PTR constant
#include "../common/transport.h"

kernel void monte_carlo(device const float *data [[buffer(0)]],
                        device atomic_float *tallies [[buffer(1)]],
                        device atomic_uint *failure [[buffer(2)]],
                        constant uint *controls [[buffer(3)]],
                        uint lane [[thread_position_in_grid]]) {
    transport_photon(data, tallies, failure, controls, lane);
}
