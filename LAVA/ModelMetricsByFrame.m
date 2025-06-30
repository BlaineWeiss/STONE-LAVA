function [SkeletonModel,StartModel,EndModel,indpointst,endpointst,r1,r2,avgr,pathlength,areasegmentsamps,rdiststores,areasections,lenareasec] = ModelMetricsByFrame(aDialine2,coordx,coordy,scaleum)
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
tic
parfor stack = 1 : size(aDialine2,3)
    %clear indarr a b rcood bdist arr midline startx
    rdistum = [];
    rdist = [];
    arr = [];
    arr2 = [];
    starty = [];
    startx = [];
    midy = [];
    midx = [];
   
    endx = [];
    endy = [];

    testbinn = aDialine2(:,:,stack);    %%%%Maybe shouldn't be an app variable (only 1 frame at a time)
    rpix = sum(testbinn,1);

    indexstart = edge(testbinn);
    indpoint = zeros(length(indexstart),1);  %Preallocation
    endpoint = zeros(length(indexstart),1);
    for i = 1:size(indexstart,2)
        for ii = 1:100
            if indexstart(ii,i) == 1
                indpoint(i) = ii;
                break
            end
        end
        for ii = 100: -1: 1
            if indexstart(ii,i) == 1
                endpoint(i) = ii;
                break
            end
        end
    end
    indpoint(1) = indpoint(2);
    endpoint(1) = endpoint(2);
    %[indpoint,indpoint2] = find(indexstart == 2)  %% indpoint is it!
    indpointst(:,stack) = indpoint;
    endpointst(:,stack) = endpoint;

    midline = floor(indpoint + (endpoint - indpoint)/2); %rpix/2);
    midline(1) = midline(2);
    cancind = find(~midline);
    midline(cancind) = [];
    indpoint(cancind) = [];
    endpoint(cancind) = [];
    for i = 1 : length(midline)
        midx(i) = floor(coordx(midline(i),i));
        midy(i) = floor(coordy(midline(i),i));

        startx(i) = round(coordx(indpoint(i),i)); %OPTION 2: ADD STACK IN 3RD DIMs
        starty(i) = round(coordy(indpoint(i),i));

        endx(i) = round(coordx(endpoint(i),i));
        endy(i) = round(coordy(endpoint(i),i));
        %OPTION 2:
        %   if coordx(indpoint(i),i) == 0
        %  startx(i) = floor(coordx(indpoint(i)+1,i,stack));
        %  starty(i) = floor(coordy(indpoint(i)+1,i,stack));
        %  end
    end

    midpairs = cat(1,midx,midy)';
    [a,ia,ic] = unique(midpairs,'rows','stable');   %%%midpoints
    areasegmentsamps{1,stack} = ic;
    startpairs = cat(1,startx,starty)';
    b = startpairs(ia,:);
    endpairs = cat(1,endx,endy)';
    be = endpairs(ia,:);
    p = pdist([b(:,1),b(:,2);be(:,1),be(:,2)]);
    x = cat(1,a(:,1),b(:,1));
    y = cat(1,a(:,2),b(:,2));
    asiz = 1:length(a);
    xpoly = fit(asiz',a(:,1),'smoothingspline','SmoothingParam',0.001);
    ypoly = fit(asiz',a(:,2),'smoothingspline','SmoothingParam',0.001);
    xs = feval(xpoly,asiz(1:1:length(asiz)));
    ys = feval(ypoly,asiz(1:1:length(asiz)));
    dist = zeros(length(ys),1);
    pnts = size(xs,1);
    totaldist = 0;
    len = 1;
    chkind = zeros(pnts, 1);
    wtc = 1;
    while wtc == 1    %Iterative algorithm for ordering skeleton model points
        chkind(len) = 1;
        if sum(chkind) == pnts
            wtc = 0;
            break
        end
        [p,I] = pdist2([xs(chkind == 0), ys(chkind ==0)], [xs(len),ys(len)] ,  'euclidean',"Smallest",1);   %distance of all set points to xs(len) ys(len
        indpnts = I;
        totaldist = totaldist + p; %pdist;
        if indpnts(1) < len - sum(chkind(1:len))  %if difference between index position and sum of points behind it is larger than the current index in chkind (some are now marked "visited").
            len = sum(chkind(1:len)) - 1 + indpnts;
        else
            ref = find(chkind == 0);
            len = ref(indpnts);
        end
    end
    for len1 = 2: length(ys)
        p = pdist([xs(len1),ys(len1) ; xs(len1-1), ys(len1-1)]);
        dist(len1-1) = p;   %%% distance between middle points in crosssections (in pixels)
    end
    dist(len1) = p;
    rdiststore3 = zeros(1,length(ys));
    rdist3 = 0;
    for len1 = 1 : length(ys) %was a

        c = [xs(len1),ys(len1) ; b(len1,1),b(len1,2)];
        ce = [xs(len1),ys(len1) ; be(len1,1),be(len1,2)];
        ar = pdist(c);     %%%Get the distance of particular midpoint a from all vessel edges.. we will then take the min to get it's closest neighbor
        ar2 = pdist(ce);
        arr(len1) = ar;
        arr2(len1) = ar2;
        %this is the fig
        %{
                        scatter(app.UIAxes,xs,ys)
                        scatter(app.UIAxes,b(:,1),b(:,2))
                        hold(app.UIAxes,'on')
                        scatter(app.UIAxes,be(:,1),be(:,2))
                        plot(app.UIAxes,c(:,1),c(:,2))
                        plot(app.UIAxes,ce(:,1),ce(:,2))
        %}
        %hold on
        %plot(app.UIAxes,ce(:,1),ce(:,2))
        %    drawnow limitrate
        %end

        rdist = arr(len1);
        rdist2 = arr2(len1);
        %rdist = bdist;          %%%% OUR RADIUS!!
        %rdist2 = bdist2;
        rdist3 = mean([rdist,rdist2]);

        %%%HERE I CAN EXTRACT THE RADIUS IN LINEAR PIXELS & USE IT TO MOVE THE ENDFOOT ROI IN THE STACK
        rdiststore3(1,len1) = rdist3;
    end

    SkeletonModel{stack} = [xs,ys];
    StartModel{stack} = b;  %coords of start and end of vessel
    EndModel{stack} = be;
    r1{stack} = arr;
    r2{stack} = arr2;
    avgr{stack} = rdist3;

    rdistum = rdiststore3 .* scaleum;
    area = rdistum .^2 * pi;
    area = filloutliers(area,"linear",'movmean',3);
    areasection = area;
    pathlength(stack) = totaldist; %sum(dist);   %%%YES!!
    rdiststores(stack,:) = {rdiststore3};    %used to approximate motion of endfoot roi
    areasections(stack,:) = {areasection};
    arlen = size(areasection,2); %was length(areasection)
    lenareasec(stack) = arlen; %number of crosssections
end
toc
end
