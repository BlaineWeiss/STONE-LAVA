function chnroi = generateROIStackImage(pixelcoords, EndfeetonChannelMask, dims)
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
datax = dims(1);
datay = dims(2);
numFrames = length(pixelcoords);
chnroi = zeros(datax, datay, numFrames, 'like', EndfeetonChannelMask);
coordlens = cellfun(@length, pixelcoords);

for i = 1:numFrames
    pix = pixelcoords{i};
    if isempty(pix)
        continue;
    end
    ind = sub2ind([datax, datay], pix(:,2), pix(:,1));
    chnroiSlice = chnroi(:,:,i);
    chnroiSlice(ind) = EndfeetonChannelMask(1:coordlens(i), i);
    chnroi(:,:,i) = chnroiSlice;
end
end
