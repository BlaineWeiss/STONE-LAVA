function [Allvalid,Valid,SkeletonModel,StartModel,EndModel,indpointst,endpointst,r1,r2,avgr,pathlength,indexingsections,areasegmentsamps,rdiststores,areasections,lenareasec] = ModelMetricsByFrame(aDialine2,coordx,coordy,CrossSectionNumber,scaleum)
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
rdistum = [];
rdist = [];
arr = [];
arr2 = [];
starty = [];
startx = [];
midy =       NaN(size(aDialine2, 2),size(aDialine2, 3));
midx =  NaN(size(aDialine2, 2),size(aDialine2, 3));
endx = [];
endy = [];





if mod(CrossSectionNumber,2) == 0
    filtnum = CrossSectionNumber + 1;
else
    filtnum = CrossSectionNumber;
end
nRows = size(coordx,1);
nCols = size(coordx,2);
Allvalid = true(nCols,1);
for stack = 1 : size(aDialine2,3)
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
        %if indpoint(i) == 0, indpoint(i) = indpoint(i-1); , end
        for ii = 100: -1: 1
            if indexstart(ii,i) == 1
                endpoint(i) = ii;
                break
            end
        end
        %if endpoint(i) == 0, endpoint(i) = endpoint(i-1); , end

    end
    indpoint(1) = indpoint(2);
    endpoint(1) = endpoint(2);
    %[indpoint,indpoint2] = find(indexstart == 2)  %% indpoint is it!
    indpointst(:,stack) = indpoint;
    endpointst(:,stack) = endpoint;

    midline = floor(indpoint + (endpoint - indpoint)/2); %rpix/2);
    midline(1) = midline(2);
    cancind = find(~midline);
    midline(cancind) = NaN; %[];
    indpoint(cancind) = NaN; %[];
    endpoint(cancind) = NaN; %[];

    %{
    for i = 1 : length(midline)
        midx(i) = floor(coordx(midline(i),i));
        midy(i) = floor(coordy(midline(i),i));

        startx(i) = round(coordx(indpoint(i),i)); %OPTION 2: ADD STACK IN 3RD DIMs
        starty(i) = round(coordy(indpoint(i),i));

        endx(i) = round(coordx(endpoint(i),i));
        endy(i) = round(coordy(endpoint(i),i));
    end
    %}
% Find valid midline entries
cols  = 1:length(midline);

valid = midline > 10 & midline <= nRows & ~isnan(midline) & indpoint > 5 & indpoint < 50 & endpoint < 95 & endpoint >50;
Allvalid = Allvalid & valid;
%midline(:) to midline(valid)
%GAP COMPATABLITY
%Unevenly distributes NOT GOOD
%{ 
    midx(valid,stack) = floor(coordx(sub2ind(size(coordx), midline(valid)', cols(valid))))';
    midy(valid,stack) = floor(coordy(sub2ind(size(coordy), midline(valid)', cols(valid))))';
    startx(valid,stack) = floor(coordx(sub2ind(size(coordx), indpoint(valid)', cols(valid))))';
    starty(valid,stack) = floor(coordy(sub2ind(size(coordy), indpoint(valid)', cols(valid))))';
    endx(valid,stack) = floor(coordx(sub2ind(size(coordx), endpoint(valid)', cols(valid))))';
    endy(valid,stack) = floor(coordy(sub2ind(size(coordy), endpoint(valid)', cols(valid))))';
%}

validmidlines = ~isnan(midline);
% Maintains Linearity of samples
%{
    midx(:,stack) = floor(coordx(sub2ind(size(coordx), midline(:)', 1:length(midline))))';
    midy(:,stack) = floor(coordy(sub2ind(size(coordy), midline(:)', 1:length(midline))))';
    startx(:,stack) = floor(coordx(sub2ind(size(coordx), indpoint(:)', 1:length(midline))))';
    starty(:,stack) = floor(coordy(sub2ind(size(coordy), indpoint(:)', 1:length(midline))))';
    endx(:,stack) = floor(coordx(sub2ind(size(coordx), endpoint(:)', 1:length(midline))))';
    endy(:,stack) = floor(coordy(sub2ind(size(coordy), endpoint(:)', 1:length(midline))))';
%}
% Can handle NaN indices, no vessel in model portions % maintains linearity of samples
    midx(validmidlines,stack) = floor(coordx(sub2ind(size(coordx), midline(validmidlines)', 1:sum(validmidlines))))';
    midy(validmidlines,stack) = floor(coordy(sub2ind(size(coordy), midline(validmidlines)', 1:sum(validmidlines))))';
    startx(validmidlines,stack) = floor(coordx(sub2ind(size(coordx), indpoint(validmidlines)', 1:sum(validmidlines))))';
    starty(validmidlines,stack) = floor(coordy(sub2ind(size(coordy), indpoint(validmidlines)', 1:sum(validmidlines))))';
    endx(validmidlines,stack) = floor(coordx(sub2ind(size(coordx), endpoint(validmidlines)', 1:sum(validmidlines))))';
    endy(validmidlines,stack) = floor(coordy(sub2ind(size(coordy), endpoint(validmidlines)', 1:sum(validmidlines))))';
%}
    %midline(:,stack) = midline;
    %indpoint(:,stack) = indpoint;
    %endpoint(:,stack) = endpoint;
end

%Allow removal of ends from volume
cutnomodelends = sum(midx,2,'includenan');
cutnomodelends = isnan(cutnomodelends);
midx = midx(~cutnomodelends,:);
midy = midy(~cutnomodelends,:);
startx = startx(~cutnomodelends,:);
starty = starty(~cutnomodelends,:);
endx = endx(~cutnomodelends,:);
endy = endy(~cutnomodelends,:);
%%%%%%%%%%%%%%%%%


for stack = 1 : size(aDialine2,3)
    midpa = cat(2,midx(:,stack),midy(:,stack));

    midpa = round(sgolayfilt(midpa, 3, filtnum));

    %DOWNSAMPLE FOR CONSISTANT CROSSSECTIONS AND EASY VOLUME
    d = [0; cumsum(sqrt(sum(diff(midpa).^2, 2)))];
    d_norm = d / d(end);
    % Remove invalid (branched regions)
    %d_norm = d_norm(valid);
    % Remove duplicates
    [d_unique, ia, ~] = unique(d_norm, 'stable');
    midpairs_unique = midpa(ia, :);
    % Define query points
    query = linspace(0, 1, CrossSectionNumber); %change to slider value%%%
    % Interpolate along the cleaned arc length
    xq = interp1(d_unique, midpairs_unique(:,1), query, 'linear');
    yq = interp1(d_unique, midpairs_unique(:,2), query, 'linear');
    resampledPoints = [xq(:), yq(:)];
    % Find closest original indices
    %idxOriginal = zeros(CrossSectionNumber, 1);
    for i = 1:CrossSectionNumber
        diffs = (midpa(:,1) - xq(i)).^2 + (midpa(:,2) - yq(i)).^2;
        [~, idxOriginal(i)] = min(diffs);
    end
    midpairs{stack} = midpa;
    aa{stack} = resampledPoints;
    idxOriginal = idxOriginal(idxOriginal > 0);
    indexingsections(:,stack) = idxOriginal; %ia;  %used for consistant cross sections (ranked analysis)
end
indexingsections = round(mean(indexingsections,2));
Valid = Allvalid(indexingsections);

areasegmentsamps = zeros(size(midpa,1), size(aDialine2,3));
parfor stack = 1 : size(aDialine2,3)
    midpairss = midpairs{stack};
    a = aa{stack};
    ia = indexingsections                           %Downsample index for model reduction and volume measurement
    [~,~,ic] = unique(midpairss,'rows','stable');   %%%midpoints
    areasegmentsamps(:,stack) = ic;                 %Preserved for max upsampling for fluorescence pairing
    startpairs = cat(2,startx(:,stack),starty(:,stack));
    %startpairs = cat(1,startx,starty)';
    b = startpairs(ia,:);
    endpairs = cat(2,endx(:,stack),endy(:,stack));
    %endpairs = cat(1,endx,endy)';
    be = endpairs(ia,:);
    p = pdist([b(:,1),b(:,2);be(:,1),be(:,2)]);
    x = cat(1,a(:,1),b(:,1));
    y = cat(1,a(:,2),b(:,2));
    asiz = 1:length(a);
    xpoly = fit(asiz',a(:,1),'smoothingspline','SmoothingParam',1);
    ypoly = fit(asiz',a(:,2),'smoothingspline','SmoothingParam',1);
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
    arr = zeros(1, length(xs));
    arr2 = zeros(1, length(xs));

    for len1 = 1 : length(ys) %was a
        c = [xs(len1),ys(len1) ; b(len1,1),b(len1,2)];
        ce = [xs(len1),ys(len1) ; be(len1,1),be(len1,2)];
        ar = pdist(c);     %%%Get the distance of particular midpoint a from all vessel edges.. we will then take the min to get it's closest neighbor
        ar2 = pdist(ce);
        arr(len1) = ar;
        arr2(len1) = ar2;
        rdist = arr(len1);
        rdist2 = arr2(len1);
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
    areasections(stack,:) = areasection;
    arlen = size(areasection,2); %was length(areasection)
    lenareasec(stack) = arlen; %number of crosssections
end
toc
end
