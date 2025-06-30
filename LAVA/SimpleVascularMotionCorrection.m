function yyy2 = SimpleVascularMotionCorrection(indpoint,endpoint,yyy,vertices2,refradvert)
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
Pixradius = floor((endpoint - indpoint)/2); %THIS IS THE RADIUS IN LINEARIZED MODEL PIXELS
midpnt = floor((indpoint + endpoint)/2);
yyy2 = repmat(yyy,1,size(endpoint,2));
vertexradius = Pixradius(vertices2(:,1),:);
midpnts = midpnt(vertices2(:,1),:);

indexfullposshift = vertices2(:,2) > endpoint(vertices2(:,1),:); %Get boundary of roi that is beyond vessel boundary
yyy2(indexfullposshift) = min(yyy2(indexfullposshift) + (vertexradius(indexfullposshift) - refradvert(indexfullposshift))./2,100);
indexfullnegshift = vertices2(:,2) < indpoint(vertices2(:,1),:);
yyy2(indexfullnegshift) = max(yyy2(indexfullnegshift) - (vertexradius(indexfullnegshift) - refradvert(indexfullnegshift))./2,1);

remaining = ~(indexfullnegshift + indexfullposshift);
yyy2(remaining) = max(min(yyy2(remaining) + (vertexradius(remaining) - refradvert(remaining)).*((yyy2(remaining)-midpnts(remaining))./vertexradius(remaining)),100),1);

end
