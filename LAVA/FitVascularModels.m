function [figfX, figfY, figf2X, figf2Y] = FitVascularModels(pos_response,posres2)
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

xx = 1:length(pos_response(:,1));
xmodel = pos_response(:,1);
%fx1coeff = polyfit(xx,xmodel,5)             %New Polyfit
[f1x, ~, output] = fit(xx',xmodel,'smoothingspline');   %%%Using FIT     CURRENTLY IN USE

yy = 1:length(pos_response(:,2));
ymodel = pos_response(:,2);
%fy1coeff = polyfit(yy,ymodel,5)             %New Polyfit
f1y = fit(yy',ymodel,'smoothingspline');  %%%Using FIT      CURRENTLY IN USE

xx2 = 1:length(posres2(:,1));
x2model = posres2(:,1);
%fx2coeff = polyfit(xx2,x2model,5)           %New Polyfit
f2x = fit(xx2',x2model,'smoothingspline')  %%%Using FIT    CURRENTLY IN USE

yy2 = 1:length(posres2(:,2));
y2model = posres2(:,2);
%fy2coeff = polyfit(yy2,y2model,5)           %New Polyfit
f2y = fit(yy2',y2model,'smoothingspline');   %%%Using FIT  CURRENTLY IN USE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIT   but with parametric equations
k=1

x1 = 1:length(pos_response(:,1));
xx = pos_response(:,1);
xv = linspace(x1(1),x1(end),600);
xx = interp1(x1,xx,xv,'pchip');

y1 = 1:length(pos_response(:,2));
yy = pos_response(:,2);
yv = linspace(y1(1),y1(end),600)
yy = interp1(y1,yy,yv,'pchip');

x2 = 1:length(posres2(:,1));
xx2 = posres2(:,1);
xv2 = linspace(x2(1),x2(end),600);
xx2 = interp1(x2,xx2,xv2,'pchip');

y2 = 1:length(posres2(:,2));
yy2 = posres2(:,2);
yv2 = linspace(y2(1),y2(end),600);
yy2 = interp1(y2,yy2,yv2,'pchip');

figfX = f1x(xv);%xx;
figfY = f1y(yv);%yy;
figf2X = f2x(xv2);%xx2;
figf2Y = f2y(yv2);% yy2;
end
