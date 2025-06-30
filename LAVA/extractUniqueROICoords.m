function outputwodups = extractUniqueROICoords(ROIcoordx, ROIcoordy)
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
[numPoints, numSweeps, numFrames] = size(ROIcoordx);
outputwodups = cell(numFrames, numSweeps);
for i = 1:numSweeps
    Sweepix = squeeze(ROIcoordx(:, i, :));
    Sweepiy = squeeze(ROIcoordy(:, i, :));

    for sta = 1:numFrames
        aa = floor([Sweepix(:, sta), Sweepiy(:, sta)]);
        aa = aa(all(aa, 2) & ~any(isnan(aa), 2), :);

        if ~isempty(aa)
            outputwodups{sta, i} = int16(unique(aa, 'rows', 'stable'));
        else
            outputwodups{sta, i} = int16([]);
        end
    end
end
end
