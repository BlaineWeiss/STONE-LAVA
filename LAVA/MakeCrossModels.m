function [slopefigf,slopefigf2,Interceptf1,Interceptf2,Invslopefigf,sysX,sysY] = MakeCrossModels(figfX,figfY,figf2X,figf2Y)
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

[sizefigfX,~] = size(figfX);
[sizefigf2X,~] = size(figf2X);

if sizefigfX > sizefigf2X
    figfX(end - (sizefigfX-sizefigf2X-1):end) = [];
    figfY(end - (sizefigfX-sizefigf2X-1):end) = [];
elseif sizefigf2X > sizefigfX
    figf2X(end - (sizefigf2X-sizefigfX-1):end) = [];
    figf2Y(end - (sizefigf2X-sizefigfX-1):end) = [];
end

dfigfX = figfX(2:length(figfX))-figfX(1:length(figfX)-1);
dfigfY = figfY(2:length(figfY))-figfY(1:length(figfY)-1);
dfigf2X = figf2X(2:length(figf2X))-figf2X(1:length(figf2X)-1);
dfigf2Y = figf2Y(2:length(figf2Y))-figf2Y(1:length(figf2Y)-1);
slopefigf2 = dfigf2Y./dfigf2X;
slopefigf = dfigfY./dfigfX;
Invslopefigf = -1./slopefigf; %perpendicular line
if size(figf2X) > size(Invslopefigf) == [1,0]
    Invslopefigf(length(Invslopefigf)+1) = Invslopefigf(length(Invslopefigf));
    slopefigf2(length(slopefigf2)+1) = slopefigf2(length(slopefigf2));
end
Interceptf1 = figfY - (figfX.*Invslopefigf); %Interceptf1 is the inverse intercept %%%need +1?
Interceptf2 = figf2Y(1:length(figf2Y)) - (figf2X(1:length(figf2X)).*slopefigf2); %%%need +1?
sysX = (Interceptf2 - Interceptf1)./(Invslopefigf - slopefigf2);   %%%Not exact.. can I fit between model and drawn fit... AND prevent overlap/cross/reuse of points along vessel
sysY = (Invslopefigf.*sysX)+Interceptf1;

end
