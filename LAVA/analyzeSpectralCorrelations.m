% Function 7: Perform spectral filtering and correlation analysis
function results = analyzeSpectralCorrelations(Parallelregion, resizedcrosssections, Times, acqrate)
    % STONE / LAVA - Scientific Analysis Software
%
% Copyright © 2025 Blaine Everett Weiss, University of Kentucky Research Foundation
% All rights reserved.
%
% This software is distributed for non-commercial academic research use only.
% Use in any commercial setting, including for-profit institutions or fee-for-service labs,
% is strictly prohibited without a separate commercial license.
%
% A patent has been filed on the underlying method. Commercial use may infringe on that IP.
% Patent application number: Pending
%
% License terms: See LICENSE_ACADEMIC.txt and LICENSE_COMMERCIAL.txt for details.
%
% For commercial licensing inquiries, contact: blaine.weiss@uky.edu
%
% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

% If you use this software in published work, please cite:

% Disrupted calcium dynamics in reactive astrocytes occur with endfeet-arteriole decoupling in an amyloid mouse model of Alzheimer’s disease
% Blaine E. Weiss, et al. bioRxiv 2025.01.24.634584; doi: https://doi.org/10.1101/2025.01.24.634584
%{
    lowcutofffreq = 0.02 / (acqrate / 2);
    highcutofffreq = [0.75 11] / (acqrate / 2);
    [b1, a1] = butter(1, lowcutofffreq, 'high');
    [b2, a2] = butter(3, [highcutofffreq(1), min(highcutofffreq(2), 0.99)]);

    numCols = desired_size;
    numFrames = length(Times);

    results.xrecsum = zeros(1, numFrames);
    results.ultralowbandstopsum = zeros(1, numFrames);
    results.VascularRegionCorr = NaN(1, numCols);
    results.VascularRegionHighfreqCorr = NaN(1, numCols);
    results.VascularRegionLowfreqCorr = NaN(1, numCols);
    results.Highbandmat = NaN(numCols, numFrames);
    results.Lowbandmat = NaN(numCols, numFrames);
    results.ROItrimmedCrossAreas = NaN(numFrames, numCols);

    for col = 1:numCols
        normarea = resizedcrosssections(:,col);
        aa = Parallelregion(col,:);
        if all(isnan(aa))
            continue
        end
        aa(isnan(aa)) = 0;

        xrec = filtfilt(b2,a2,aa);
        ultralowbandstop = filtfilt(b1,a1,aa);

        xrecEnv = envelope(xrec,3,'analytic');
        ultralowEnv = envelope(ultralowbandstop,3,'analytic');

        results.xrecsum = results.xrecsum + xrec;
        results.ultralowbandstopsum = results.ultralowbandstopsum + ultralowbandstop;
        results.VascularRegionCorr(col) = corr2(normarea, ultralowbandstop');
        results.VascularRegionHighfreqCorr(col) = corr(normarea, xrecEnv');
        results.VascularRegionLowfreqCorr(col) = corr(normarea, ultralowEnv');

        results.ROItrimmedCrossAreas(:,col) = normarea;
        results.Highbandmat(col,:) = xrecEnv;
        results.Lowbandmat(col,:) = ultralowEnv;
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
profile on;
lowcutofffreq = 0.02 / (acqrate / 2);
highcutofffreq = [0.75 11] / (acqrate / 2);
[b1, a1] = butter(1, lowcutofffreq, 'high');
[b2, a2] = butter(3, [highcutofffreq(1), min(highcutofffreq(2), 0.99)]);
numCols = size(resizedcrosssections,2);
numFrames = length(Times);

% Preallocate with NaNs for clarity and memory
results = struct();
results.ButterCoeff = NaN(1,4);
results.xrecsum = zeros(1, numFrames);
results.ultralowbandstopsum = zeros(1, numFrames);
results.VascularRegionCorr = NaN(1, numCols);
results.VascularRegionHighfreqCorr = NaN(1, numCols);
results.VascularRegionLowfreqCorr = NaN(1, numCols);
results.Highbandmat = NaN(numCols, numFrames);
results.Lowbandmat = NaN(numCols, numFrames);
results.ROItrimmedCrossAreas = NaN(numFrames, numCols);
ButterCoeff = [b1,a1,b2,a2];
NaNindex = zeros(1,numCols);

parfor col = 1:numCols
    normarea = resizedcrosssections(:,col);
    aa = Parallelregion(col,:);

    if all(isnan(aa))
        NaNindex(col) = 1;
        continue
    end

    aa(isnan(aa)) = 0;

    try
        xrec = filtfilt(b2, a2, aa);
        ultralowbandstop = filtfilt(b1, a1, aa);
    catch
        continue;
    end

    xrecEnv = envelope(xrec,3,'analytic');
    ultralowEnv = envelope(ultralowbandstop,3,'analytic');

    xrecsum_temp(col,:) = xrec;
    ultralowbandstopsum_temp(col,:) = ultralowbandstop;
    VascularRegionCorr(col) = corr(normarea, ultralowbandstop', 'rows','complete');
    VascularRegionHighfreqCorr(col) = corr(normarea, xrecEnv', 'rows','complete');
    VascularRegionLowfreqCorr(col) = corr(normarea, ultralowEnv', 'rows','complete');

    ROItrimmedCrossAreas(:,col) = normarea;
    Highbandmat(col,:) = xrecEnv;
    Lowbandmat(col,:) = ultralowEnv;
end
xrecsum = sum(xrecsum_temp, 1);
ultralowbandstopsum = sum(ultralowbandstopsum_temp, 1);

% Package the results
results = struct(...
    'ButterCoeff', ButterCoeff, ...
    'xrecsum', xrecsum, ...
    'ultralowbandstopsum', ultralowbandstopsum, ...
    'VascularRegionCorr', VascularRegionCorr, ...
    'VascularRegionHighfreqCorr', VascularRegionHighfreqCorr, ...
    'VascularRegionLowfreqCorr', VascularRegionLowfreqCorr, ...
    'Highbandmat', Highbandmat, ...
    'Lowbandmat', Lowbandmat, ...
    'ROItrimmedCrossAreas', ROItrimmedCrossAreas, ...
    'NaNindex', NaNindex);
profile viewer
end
