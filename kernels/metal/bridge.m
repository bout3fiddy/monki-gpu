// Local Metal host bridge. Physics is exclusively in transport.metal.
// ARC owns Metal resources through one opaque session retained by Python.
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include <stdint.h>
#include <string.h>

@interface MonkiMetalSession : NSObject
@property id<MTLDevice> device;
@property id<MTLCommandQueue> queue;
@property id<MTLComputePipelineState> pipeline;
@property id<MTLBuffer> scene;
@property id<MTLBuffer> tallies;
@property id<MTLBuffer> failure;
@property double gpuSeconds;
@property size_t roundBytes;
@end
@implementation MonkiMetalSession
@end

static void error_text(char *buffer, size_t capacity, NSString *message) {
    snprintf(buffer, capacity, "%s", (message ?: @"Unknown Metal error").UTF8String);
}

// Version the ctypes ABI so stale bridges fail before any mismatched native call.
int monki_metal_abi_version(void) {
    return 2;
}

void *monki_metal_create(const char *source, uint32_t device_index, const float *scene,
                         size_t scene_count, size_t tally_count, uint32_t round_capacity,
                         char *error, size_t error_size) {
    @autoreleasepool {
        NSArray<id<MTLDevice>> *devices = MTLCopyAllDevices();
        if (device_index >= devices.count) {
            error_text(error, error_size, @"Requested Metal device is unavailable");
            return NULL;
        }

        MonkiMetalSession *session = [MonkiMetalSession new];
        session.device = devices[device_index];
        if (![session.device supportsFamily:MTLGPUFamilyMetal3]) {
            error_text(error, error_size, @"Metal 3 floating-point atomics are required");
            return NULL;
        }
        session.queue = [session.device newCommandQueue];
        MTLCompileOptions *options = [MTLCompileOptions new];
        options.fastMathEnabled = NO;
        options.languageVersion = MTLLanguageVersion3_0;
        NSError *failure = nil;
        id<MTLLibrary> library =
            [session.device newLibraryWithSource:[NSString stringWithUTF8String:source]
                                         options:options
                                           error:&failure];
        if (!library) {
            error_text(error, error_size, failure.description);
            return NULL;
        }

        id<MTLFunction> transport = [library newFunctionWithName:@"monte_carlo"];
        if (!transport) {
            error_text(error, error_size, @"Metal kernel is missing monte_carlo");
            return NULL;
        }
        session.pipeline = [session.device newComputePipelineStateWithFunction:transport
                                                                         error:&failure];
        session.scene = [session.device newBufferWithBytes:scene
                                                    length:scene_count * sizeof(float)
                                                   options:MTLResourceStorageModeShared];
        session.roundBytes = tally_count * sizeof(float);
        session.tallies = [session.device newBufferWithLength:session.roundBytes * round_capacity
                                                      options:MTLResourceStorageModeShared];
        session.failure = [session.device newBufferWithLength:sizeof(uint32_t)
                                                      options:MTLResourceStorageModeShared];
        if (!session.queue || !session.pipeline || !session.scene || !session.tallies ||
            !session.failure) {
            error_text(error, error_size,
                       failure ? failure.description : @"Metal allocation failed");
            return NULL;
        }

        return (__bridge_retained void *)session;
    }
}

// One bounded dispatch, copied out only after GPU completion. A nonzero shader
// failure invalidates the round; never silently accept truncated photon histories.
int monki_metal_run(void *handle, const uint32_t *controls, float *output, char *error,
                    size_t error_size) {
    @autoreleasepool {
        MonkiMetalSession *session = (__bridge MonkiMetalSession *)handle;
        size_t activeBytes = session.roundBytes * controls[7];
        if (!controls[6] || !controls[7] || activeBytes > session.tallies.length) {
            error_text(error, error_size, @"Invalid Metal round group");
            return 1;
        }
        memset(session.tallies.contents, 0, activeBytes);
        memset(session.failure.contents, 0, sizeof(uint32_t));
        id<MTLCommandBuffer> command = [session.queue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [command computeCommandEncoder];
        if (!command || !encoder) {
            error_text(error, error_size, @"Metal command allocation failed");
            return 1;
        }
        [encoder setComputePipelineState:session.pipeline];
        [encoder setBuffer:session.scene offset:0 atIndex:0];
        [encoder setBuffer:session.tallies offset:0 atIndex:1];
        [encoder setBuffer:session.failure offset:0 atIndex:2];
        [encoder setBytes:controls length:8 * sizeof(uint32_t) atIndex:3];
        NSUInteger width = MIN(64, session.pipeline.maxTotalThreadsPerThreadgroup);
        [encoder dispatchThreads:MTLSizeMake(controls[0], 1, 1)
            threadsPerThreadgroup:MTLSizeMake(width, 1, 1)];
        [encoder endEncoding];
        [command commit];
        [command waitUntilCompleted];

        if (command.status != MTLCommandBufferStatusCompleted) {
            error_text(error, error_size, command.error.description);
            return 1;
        }
        // Read command timestamps only after completion; no extra GPU synchronization.
        if (controls[5])
            session.gpuSeconds += command.GPUEndTime - command.GPUStartTime;
        uint32_t flag = *(uint32_t *)session.failure.contents;
        if (flag) {
            snprintf(error, error_size, "Metal transport failed (flag %u); round rejected", flag);
            return 1;
        }

        memcpy(output, session.tallies.contents, activeBytes);
        return 0;
    }
}

double monki_metal_gpu_seconds(void *handle) {
    MonkiMetalSession *session = (__bridge MonkiMetalSession *)handle;
    return session.gpuSeconds;
}

const char *monki_metal_device_name(void *handle) {
    MonkiMetalSession *session = (__bridge MonkiMetalSession *)handle;
    return session.device.name.UTF8String;
}

void monki_metal_close(void *handle) {
    @autoreleasepool {
        MonkiMetalSession *session = (__bridge_transfer MonkiMetalSession *)handle;
        (void)session;
    }
}
