// CUDA driver/NVRTC host adapter. Each session owns a context and resident buffers.
#include <cuda.h>
#include <nvrtc.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#ifndef MONKI_CUDA_INCLUDE
#define MONKI_CUDA_INCLUDE "/usr/local/cuda/include"
#endif

static void checked(CUresult status) {
    if (status != CUDA_SUCCESS) {
        const char *message = nullptr;
        cuGetErrorString(status, &message);
        throw std::runtime_error(message ? message : "CUDA driver failure");
    }
}
static void checked(nvrtcResult status) {
    if (status != NVRTC_SUCCESS)
        throw std::runtime_error(nvrtcGetErrorString(status));
}
static void error_text(char *error, size_t capacity, const char *message) {
    if (error && capacity)
        std::snprintf(error, capacity, "%s", message);
}
struct Session {
    CUcontext context = nullptr;
    CUmodule module = nullptr;
    CUfunction kernel = nullptr;
    CUdeviceptr scene = 0, tallies = 0, failure = 0, controls = 0;
    CUevent begin = nullptr, end = nullptr;
    size_t round_bytes = 0;
    uint32_t capacity = 0;
    double seconds = 0;
    char name[256] = {};
    char runtime_info[512] = {};
    ~Session() {
        // Context destruction also releases allocations after a partial create.
        if (context)
            cuCtxDestroy(context);
    }
};
struct ContextScope {
    explicit ContextScope(CUcontext context) {
        checked(cuCtxPushCurrent(context));
    }
    ~ContextScope() {
        CUcontext previous;
        cuCtxPopCurrent(&previous);
    }
};
struct Program {
    nvrtcProgram value = nullptr;
    ~Program() {
        if (value)
            nvrtcDestroyProgram(&value);
    }
};
extern "C" int monki_cuda_abi_version() {
    return 2;
}
extern "C" void *monki_cuda_create(const char *source, uint32_t device_index, const float *scene,
                                   size_t scene_count, size_t tally_count, uint32_t capacity,
                                   char *error, size_t error_size) {
    try {
        checked(cuInit(0));
        int count = 0;
        checked(cuDeviceGetCount(&count));
        if (device_index >= static_cast<uint32_t>(count))
            throw std::runtime_error("Requested CUDA device is unavailable");
        if (!source || !scene || !scene_count || !tally_count || !capacity || capacity > 32)
            throw std::runtime_error("Invalid CUDA session dimensions");
        auto session = std::make_unique<Session>();
        CUdevice device;
        checked(cuDeviceGet(&device, static_cast<int>(device_index)));
        checked(cuDeviceGetName(session->name, sizeof(session->name), device));
        int major, minor;
        checked(cuDeviceGetAttribute(&major, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MAJOR, device));
        checked(cuDeviceGetAttribute(&minor, CU_DEVICE_ATTRIBUTE_COMPUTE_CAPABILITY_MINOR, device));
        int compiler_major, compiler_minor, driver_version;
        checked(nvrtcVersion(&compiler_major, &compiler_minor));
        checked(cuDriverGetVersion(&driver_version));
        std::snprintf(session->runtime_info, sizeof(session->runtime_info),
                      "NVRTC %d.%d; driver %d; compute_%d%d; c++17; fmad=false; prec-div=true; "
                      "prec-sqrt=true; ftz=false",
                      compiler_major, compiler_minor, driver_version, major, minor);
        checked(cuCtxCreate(&session->context, 0, device));
        CUcontext created;
        checked(cuCtxPopCurrent(&created));
        ContextScope scope(session->context);
        Program program;
        checked(nvrtcCreateProgram(&program.value, source, "monki.cu", 0, nullptr, nullptr));
        // Keep division/sqrt precise and preserve explicit operation boundaries.
        std::string architecture =
            "--gpu-architecture=compute_" + std::to_string(major) + std::to_string(minor);
        const char *options[] = {"--std=c++17",       "--fmad=false",
                                 "--prec-div=true",   "--prec-sqrt=true",
                                 "--ftz=false",       "--include-path=" MONKI_CUDA_INCLUDE,
                                 architecture.c_str()};
        nvrtcResult compiled =
            nvrtcCompileProgram(program.value, sizeof(options) / sizeof(*options), options);
        if (compiled != NVRTC_SUCCESS) {
            size_t size;
            checked(nvrtcGetProgramLogSize(program.value, &size));
            std::vector<char> log(size);
            checked(nvrtcGetProgramLog(program.value, log.data()));
            throw std::runtime_error(std::string("CUDA compilation failed: ") + log.data());
        }
        size_t size;
        checked(nvrtcGetPTXSize(program.value, &size));
        std::vector<char> ptx(size);
        checked(nvrtcGetPTX(program.value, ptx.data()));
        checked(cuModuleLoadData(&session->module, ptx.data()));
        checked(cuModuleGetFunction(&session->kernel, session->module, "monte_carlo"));
        session->round_bytes = tally_count * sizeof(float);
        session->capacity = capacity;
        checked(cuMemAlloc(&session->scene, scene_count * sizeof(float)));
        checked(cuMemcpyHtoD(session->scene, scene, scene_count * sizeof(float)));
        checked(cuMemAlloc(&session->tallies, session->round_bytes * capacity));
        checked(cuMemAlloc(&session->failure, sizeof(uint32_t)));
        checked(cuMemAlloc(&session->controls, 8 * sizeof(uint32_t)));
        checked(cuEventCreate(&session->begin, CU_EVENT_DEFAULT));
        checked(cuEventCreate(&session->end, CU_EVENT_DEFAULT));
        return session.release();
    } catch (const std::exception &failure) {
        error_text(error, error_size, failure.what());
        return nullptr;
    }
}
extern "C" int monki_cuda_run(void *handle, const uint32_t *controls, float *output, char *error,
                              size_t error_size) {
    try {
        auto *session = static_cast<Session *>(handle);
        if (!session || !controls || !output || !controls[0] || controls[0] > 65536 ||
            !controls[6] || !controls[7] || controls[7] > session->capacity ||
            uint64_t(controls[1]) + controls[0] > uint64_t(controls[6]) * controls[7])
            throw std::runtime_error("Invalid CUDA round group");
        ContextScope scope(session->context);
        size_t active = session->round_bytes * controls[7];
        checked(cuMemsetD8(session->tallies, 0, active));
        checked(cuMemsetD8(session->failure, 0, sizeof(uint32_t)));
        checked(cuMemcpyHtoD(session->controls, controls, 8 * sizeof(uint32_t)));
        void *arguments[] = {&session->scene, &session->tallies, &session->failure,
                             &session->controls};
        if (controls[5])
            checked(cuEventRecord(session->begin, nullptr));
        checked(cuLaunchKernel(session->kernel, (controls[0] + 127) / 128, 1, 1, 128, 1, 1, 0,
                               nullptr, arguments, nullptr));
        if (controls[5])
            checked(cuEventRecord(session->end, nullptr));
        checked(cuCtxSynchronize());
        if (controls[5]) {
            float milliseconds;
            checked(cuEventElapsedTime(&milliseconds, session->begin, session->end));
            session->seconds += milliseconds / 1000.0;
        }
        uint32_t flag;
        checked(cuMemcpyDtoH(&flag, session->failure, sizeof(flag)));
        if (flag)
            throw std::runtime_error("CUDA transport failed (flag " + std::to_string(flag) +
                                     "); round rejected");
        checked(cuMemcpyDtoH(output, session->tallies, active));
        return 0;
    } catch (const std::exception &failure) {
        error_text(error, error_size, failure.what());
        return 1;
    }
}
extern "C" double monki_cuda_gpu_seconds(void *handle) {
    return static_cast<Session *>(handle)->seconds;
}
extern "C" const char *monki_cuda_device_name(void *handle) {
    return static_cast<Session *>(handle)->name;
}
extern "C" void monki_cuda_close(void *handle) {
    delete static_cast<Session *>(handle);
}

extern "C" const char *monki_cuda_runtime_info(void *handle) {
    return static_cast<Session *>(handle)->runtime_info;
}
