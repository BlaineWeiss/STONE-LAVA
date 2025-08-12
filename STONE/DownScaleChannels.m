function [] = DownScaleChannels(Arraysize,s,blocksize)
UmPix = SharedData.getData("UmPix");
Chn1 = SharedData.getData("Chn1");
Chn2 = SharedData.getData("Chn2");
FOV = SharedData.getData("FOV");
%%%This gets more complicated in the case of rectangular rois.. The pixels are still square, but the amount of blocks, and the remainders can be different.
%Need to separate rema into a remaX and a remaY
remaX = rem(Arraysize(1), s)
remaY = rem(Arraysize(2), s)
app.rema = [remaX,remaY]; %Used if I want to remap
assignin('base','rema',app.rema)
if sum(app.rema) ~= 0
    remaXSt = floor(remaX/2);
    remaXend = remaX - remaXSt;
    remaYSt = floor(remaY/2);
    remaYend = remaY - remaYSt;

    trimchn = Chn1(1+remaXSt:end-remaXend,1+remaYSt:end-remaYend,:); %trims equal excess on all sides
    trimchn2 = Chn2(1+remaXSt:end-remaXend,1+remaYSt:end-remaYend,:);
else
    trimchn = Chn1;
    trimchn2 = Chn2;
end
FOV = size(trimchn,[1 2]).* UmPix;  %app.
SharedData.setData("FOV",FOV)
pixarea = (s.^2) .* (UmPix.*UmPix); %app.
SharedData.setData("pixarea",pixarea);
T = Arraysize(3);


arraymem = whos('trimchn')
gpumemory = GPUDeviceCheck;
if gpumemory > 1 
% Upload to GPU
Chn1_gpu = gpuArray(trimchn);
Chn2_gpu = gpuArray(trimchn);
end
% Get trimmed sizes
[H_trim, W_trim, ~] = size(Chn1_gpu);
Ny = ceil(H_trim / s);
Nx = ceil(W_trim / s);

% Reshape to 5D: [blockRows, Ny, blockCols, Nx, T]
Chn1_blocks = reshape(Chn1_gpu, s, Ny, s, Nx, T);
Chn2_blocks = reshape(Chn2_gpu, s, Ny, s, Nx, T);

% Permute to [blockRows, blockCols, Ny, Nx, T]
Chn1_blocks = permute(Chn1_blocks, [1, 3, 2, 4, 5]);
Chn2_blocks = permute(Chn2_blocks, [1, 3, 2, 4, 5]);

% Compute mean over rows and cols
Chn1_ds_gpu = mean(mean(Chn1_blocks, 1), 2);  % [1, 1, Ny, Nx, T]
Chn2_ds_gpu = mean(mean(Chn2_blocks, 1), 2);

% Reshape to [Ny, Nx, T] and move back to CPU
Chn1_ds = reshape(gather(Chn1_ds_gpu), Ny, Nx, T);
Chn2_ds = reshape(gather(Chn2_ds_gpu), Ny, Nx, T);


%    fun = @(block_struct) mean(block_struct.data,[1,2])
%    Chn1_ds = blockproc(trimchn,blocksize,fun,'UseParallel',true);  %%% average pixels into blocks (sliding would be better.. would use center point for indexing) <-NEXT UPDATE!
%    Chn2_ds = blockproc(trimchn2,blocksize,fun,'UseParallel',true);
%    Chn1_ds = single(Chn1_ds);
%    Chn2_ds = single(Chn2_ds);

SharedData.setData("Chn1_ds",Chn1_ds)
SharedData.setData("Chn2_ds",Chn2_ds)
end