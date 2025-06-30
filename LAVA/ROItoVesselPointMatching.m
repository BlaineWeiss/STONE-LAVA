function [vertices2,refradvert,xxx,yyy] = ROItoVesselPointMatching(ROI,refrad)
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

% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

% If you use this software in published work, please cite:

% Disrupted calcium dynamics in reactive astrocytes occur with endfeet-arteriole decoupling in an amyloid mouse model of Alzheimer’s disease
% Blaine E. Weiss, et al. bioRxiv 2025.01.24.634584; doi: https://doi.org/10.1101/2025.01.24.634584
%%
vertices = round(ROI.Position);
numpointsvert = length(vertices);
xxx = interp(vertices(:,1),2);
yyy = interp(vertices(:,2),2);

vertices2 = [xxx,yyy];
vertices2 = round(vertices2);
vertices2(vertices2 < 1) = 1;
yyy = vertices2(:,2);

vertices2(vertices2 > length(refrad)) = length(refrad);  %points out left and right past linar model will be constrained to boundary
minvertpnt = min(vertices2(:,1));
maxvertpnt = min(max(vertices2(:,1)),length(refrad));
try
refradvert = refrad(vertices2(:,1)); %Get radius reference from vertices of roi polygon
catch
refradvert = nan(1,maxvertpnt);
refradvert(minvertpnt:maxvertpnt) = refrad(vertices2(:,1));
end

end
