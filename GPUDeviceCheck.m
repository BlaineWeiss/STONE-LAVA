function gpumemory = GPUDeviceCheck()
%reset(gpuDevice)
gpumemory = gpuDevice().AvailableMemory

end