function [pixelcoords, bound,ROIcoordx,ROIcoordy] = processROIcoords(ratOutVes, xxx, yyy2, coordx, coordy, dims)
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
for stack = 1:dims(3)
endfoot = roipoly(ratOutVes,xxx,yyy2(:,stack));
ROIcoordx(:,:,stack) = uint16(round(endfoot .* coordx));    %%%The coordinates of the ROI on the channel image
ROIcoordy(:,:,stack) = uint16(round(endfoot .* coordy));
end

%%%%%%%%%%%%%%%
wkROIcoordx = reshape(ROIcoordx,[],dims(3));
wkROIcoordy = reshape(ROIcoordy,[],dims(3));
coordinatepairs = permute(cat(3,wkROIcoordx,wkROIcoordy),[1,3,2]);
parfor stack = 1:dims(3)
[C,ia,ic] = unique(coordinatepairs(:,:,stack),'rows');  %%%% I think that is it!!! size(C,1) is the amount of pixels and thus the REAL pixel area (ia = input index & ic= output index)
Thepixels = coordinatepairs(ia,:,stack);  
%%%ELIMINATE 0,0 COORDINATE
Thepixels(1,:) = [];

[rex,~] = find(Thepixels(:,2) > dims(1));
Thepixels(rex,:) = [];   %nulls out of image bounds pixels
[rex,~] = find(Thepixels(:,1) > dims(2));
Thepixels(rex,:) = [];   %nulls out of image bounds pixels


if size(isnan(Thepixels),1) > 0 
   excludeam = isnan(Thepixels);
   excludeam = sum(excludeam,2);
   Thepixels(excludeam > 0,:) = [];
end
pixelcoords1{stack,1} = double(Thepixels);  %Added 5/25/23 to conserve roi in motion
pixelcoords{stack,1} = single(Thepixels);
bound{stack} = uint16(floor(boundary(pixelcoords1{stack,1}(:,1),pixelcoords1{stack,1}(:,2),1)));
end
end
