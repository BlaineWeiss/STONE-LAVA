function EndfeetonChannelMask = extractFluorescenceTraceFromCoords(pixelcoords, Chn1)
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
coordlens = cellfun(@length, pixelcoords);
coordlen = max(coordlens);
numFrames = size(Chn1, 3);
EndfeetonChannelMask = zeros(coordlen, numFrames, 'like', Chn1);

for stack = 1:numFrames
    pix = pixelcoords{stack};
    if isempty(pix)
        continue;
    end
    inds = sub2ind([size(Chn1,1), size(Chn1,2)], pix(1:coordlens(stack),2), pix(1:coordlens(stack),1));
    EndfeetonChannelMask(1:coordlens(stack), stack) = Chn1(inds + (stack-1)*size(Chn1,1)*size(Chn1,2));
end

%{
coordlens = cellfun(@length, pixelcoords);
coordlen = max(coordlens);
numFrames = size(Chn1, 3);
EndfeetonChannelMask = zeros(coordlen, numFrames);

for stack = 1:numFrames
    pix = pixelcoords{stack};
    for k = 1:coordlens(stack)
        EndfeetonChannelMask(k, stack) = Chn1(pix(k,2), pix(k,1), stack);
    end
end
%}
end
