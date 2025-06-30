function Parallelregion = extractPerivascularIntensityBands(LinearizedCaExtract, LinearizedCaMask, ROIcoordx, ROIcoordy, icarea, parind)
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

%{
profile on;
numFrames = size(LinearizedCaExtract, 3);
Parallelregion = zeros(parind, numFrames, 'like', LinearizedCaExtract);

parfor stack = 1:numFrames
    x = ROIcoordx(:,:,stack);
    y = ROIcoordy(:,:,stack);
    Val = LinearizedCaExtract(:,:,stack);

    for col = 1:parind
        idx = find(col == icarea);
        if isempty(idx)
            continue;
        end
        xx = x(:,idx);
        yy = y(:,idx);
        Ca = Val(:,idx);

        xx1 = reshape(xx,[],1);
        yy1 = reshape(yy,[],1);
        Ca1 = reshape(Ca,[],1);

        coords = cat(2,yy1,xx1,Ca1);
        coordsext = unique(coords, 'rows', 'stable');
        if ~isempty(coordsext)
            coordsext(1,:) = []; % Remove (0,0) or duplicate pixel artifact
            if ~isempty(coordsext)
                Parallelregion(col,stack) = mean(coordsext(:,3), 'omitnan');
            else
                Parallelregion(col,stack) = NaN;
            end
        else
            Parallelregion(col,stack) = NaN;
        end
    end
end
profile viewer
%}

profile on;
numFrames = size(LinearizedCaExtract, 3);
Parallelregion = NaN(parind, numFrames, 'like', LinearizedCaExtract);

% Precompute icarea mapping to reduce repeated calls in loop
map = cell(parind, 1);
for col = 1:parind
    map{col} = find(icarea == col);
end

parfor stack = 1:numFrames
    x = ROIcoordx(:,:,stack);
    y = ROIcoordy(:,:,stack);
    Ca = LinearizedCaExtract(:,:,stack);

    for col = 1:parind
        idx = map{col};
        if isempty(idx)
            continue;
        end

        xx = x(:,idx);
        yy = y(:,idx);
        C  = Ca(:,idx);

        coords = [reshape(yy,[],1), reshape(xx,[],1), reshape(C,[],1)];
        if isempty(coords)
            continue;
        end

        coordsext = unique(coords, 'rows', 'stable');
        if size(coordsext, 1) > 1
            Parallelregion(col,stack) = mean(coordsext(2:end,3), 'omitnan');
        end
    end
end
profile viewer
end
