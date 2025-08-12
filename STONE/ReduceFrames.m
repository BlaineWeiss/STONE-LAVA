function averagedFrameschn = ReduceFrames(Chn,input,Duration) %Chn1 - Chn ; in - input  ; seconds - Duration

numframes = size(Chn,3);
secondstoplaythrough = numframes./ (in(1) .* in(2));  %how long will video take to playthrough
startingframeres = numframes ./ Duration;   %In case the acq frameaveraged
rollingframeavgf = max(ceil(startingframeres ./ in(3)),1);  %without thresholding the lower end to 1 im(3) if beyond stack resolution, will upsample the stack
Xplayback = startingframeres.* in(1);
ratioadjustmentrate = Xplayback ./ in(2) ;
ratioadjustmentrate = max(ratioadjustmentrate,1); %if less than 1, then it will not downscale but instead upscale(duplicateframes)
smoothedchn = movmean(Chn1,rollingframeavgf,3);  %was im(3)
averagedFrameschn = [];
% Loop through the frames and average every nth frame
for i = 1:ratioadjustmentrate:numframes   %increment was in(3) replaced with ratio for Xplayback buffer refresh rates
    framed = min(round(i),numframes);
    % Extract frames to be averaged
    FramesToAverage1 = smoothedchn(:, :, framed:min(framed+rollingframeavgf-1, numframes)); %was in(3) before set to sec res
    % Calculate the mean of the frames along the third dimension (time)
    averagedFrame1 = mean(FramesToAverage1, 3);
    % Append the averaged frame to the result
    averagedFrameschn = cat(3, averagedFrameschn, averagedFrame1);
end
end