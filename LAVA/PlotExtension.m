function [Y,x,Extension,xi,yi] = PlotExtension(Chn2,axisobj,UmPix,figfX,figfY,sysX,sysY,Invslopefigf,Interceptf1,Extension,passthrough)
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
g = canUseGPU
if g == 1
xi = gpuArray(zeros(1,100));  %or size(Y,2); %%%Consider switching if need higher resolution in the future (specify in x linspace function)
yi = gpuArray(zeros(1,100));
else
xi = zeros(1,100);  %or size(Y,2); %%%Consider switching if need higher resolution in the future (specify in x linspace function)
yi = zeros(1,100);
end
%app.ModelExtension = Extension
slidvalue = Extension;
crosschange = round(slidvalue ./ UmPix(1));
imagesc(axisobj,Chn2(:,:,1));
axis(axisobj);
hold(axisobj,'on');
yinv = 0;
for fi = 1:length(figfX)
    figs = figfX(fi); %%%%%%SPEED
    sysXx = sysX(fi);
    sysYy = sysY(fi);
    figys = figfY(fi);
    if figs < sysXx
        x1 = figs;
        x2 = sysXx;
        y1 = figys;
        y2 = sysYy;
    else
        x1 = sysXx;
        x2 = figs;
        y1 = sysYy;
        y2 = figys;
    end
    x = linspace(x1,x2); %Creates X range of 100 points bisect vessel**** SHIFT VALUE EXPANDS RANGE (SMALLER THE MORE REPETITIVE, TOO LARGE SKIPS PIXELS) Problem arises if roi is on edge of image... would need smaller spread
    Y = (Invslopefigf(fi).*(x))+Interceptf1(fi); %Calculates Coordinate for bisected line
    seekerx = find(x > size(Chn2,1) | x < 1);
    seekery = find(Y > size(Chn2,2)| Y < 1 );
    theta = atan(Invslopefigf(fi));
    Ymovement = abs(crosschange*sin(theta));
    if sysYy < figys
        y1 = sysYy;
        y2 = figys;
        Y = linspace(y1-Ymovement,y2+Ymovement);
    else
        y1 = figys;
        y2 = sysYy;
        Y = linspace(y1-Ymovement,y2+Ymovement);
        Y = flip(Y);
    end
    %  Y = linspace(y1-Ymovement,y2+Ymovement);  % New Y range
    x = (Y - Interceptf1(fi))./Invslopefigf(fi); %Plug Y into cross section formula.

    % Y = (Invslopefigf(fi).*(x)) + Interceptf1(fi);
    plot(x,Y,'Parent',axisobj);
    drawnow limitrate;
    switch passthrough
        case 0
            if fi ~= 1
                [xsiz,~] = size(xi);
                xi(xsiz+1,:) = single(floor(x)) ;  %X val of lines down rows %Stores X bisect line val 1
                yi(xsiz+1,:) = single(floor(Y)) ;  %Y val of lines down rows %Stores Y b isect calculation val 1
                %  NDiam(length(NDiam)+1) = plot(x,Y); %plots of lines next in iteration
            else
                xsiz = 1 ;
                xi(xsiz,:) = single(floor(x));  %Stores remaining values
                yi(xsiz,:) = single(floor(Y));
            end
        case 1
    end
end
drawnow;
end