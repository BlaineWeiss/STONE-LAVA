function [outputpixvals, averageint] = sampleAndAverage(Chn, outputwodups)
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
%%
[numFrames, numSweeps] = size(outputwodups);
outputpixvals = cell(numFrames, numSweeps);
averageint = NaN(numFrames, numSweeps);
frameSize = size(Chn{1});

parfor i = 1:numSweeps
    localPixVals = cell(numFrames, 1);
    localAvg = NaN(numFrames, 1);

    for ii = 1:numFrames
        coords = outputwodups{ii, i};
        if ~isempty(coords) && size(coords, 2) == 2
            ind = sub2ind(frameSize, coords(:, 2), coords(:, 1));
            vals = Chn{ii}(ind);

            localPixVals{ii} = vals;
            localAvg(ii) = mean(vals, 'omitnan');
        else
            localPixVals{ii} = NaN;
        end
    end

    outputpixvals(:, i) = localPixVals;
    averageint(:, i) = localAvg;
end
%{
profile on
pardowntime = size(outputwodups,1);
for i = 1:size(outputwodups, 2) % Sweep across vessels
    for ii = 1:pardowntime % Down time
        crosspixcoords = outputwodups{ii, i}; % Get coordinates to sample then average
        siz = size(crosspixcoords, 1);
        siz2 = size(crosspixcoords,2);
        if siz > 0 & siz2 > 1 % Only if there are coordinates
            a = sub2ind(size(Chn{1},[1,2]), crosspixcoords(:, 2), crosspixcoords(:, 1));
            outputpixvals{ii, i} = Chn{ii}(a);
            averageint(ii, i) = mean(a, 'omitnan'); % Average, omitting NaNs
        else
            % Handle the case where there are no coordinates
            outputpixvals{ii, i} = NaN; % Fill with NaNs if no coordinates
            averageint(ii, i) = NaN; % Set average to NaN
        end
    end
end
profile viewer
%}
end
