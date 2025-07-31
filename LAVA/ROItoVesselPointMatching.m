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
% For commercial licensing inquiries, contact: [Your Email] or [Tech Transfer Office Email]
%
% THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

% If you use this software in published work, please cite:
% [Your Citation Here]
%%
vertices = round(ROI.Position);
numpointsvert = length(vertices);

%xxx = interp(vertices(:,1),2);
%yyy = interp(vertices(:,2),2);


vertices = ROI.Position;  % Nx2 matrix: [x, y]
vertices(end+1,:) = vertices(1,:);  % Close the polygon

% Compute cumulative distance along the polygon boundary
dists = [0; cumsum(sqrt(sum(diff(vertices).^2, 2)))];

% Total boundary length and how many points you want
totalLength = dists(end);
numInterpPoints = round(size(vertices,1) * 3);  % e.g. triple the number of points

% New equally spaced distance values
distsInterp = linspace(0, totalLength, numInterpPoints);

% Interpolate X and Y separately using distance as param
xxx = interp1(dists, vertices(:,1), distsInterp, 'linear');
yyy = interp1(dists, vertices(:,2), distsInterp, 'linear');


vertices2 = [xxx',yyy'];
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