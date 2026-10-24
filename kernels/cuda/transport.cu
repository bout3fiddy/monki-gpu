// NVIDIA entry point; one CUDA thread owns one complete photon history.
#include "platform.cuh"
#include "../common/transport.h"
extern "C" __global__ void monte_carlo(const float *data, float *tallies, unsigned int *failure,
                                       const unsigned int *controls) {
    transport_photon(data, tallies, failure, controls, blockIdx.x * blockDim.x + threadIdx.x);
}
