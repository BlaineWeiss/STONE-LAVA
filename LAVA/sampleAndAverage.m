function [packedModelstructure, averageTS, LegacyCellVals] = sampleAndAverage(Chn, packedModelstructure,varargin)
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
% Inputs
%   Chn          : {T×1} cell, each Chn{t} is a [H×W] numeric frame
%   coordPacked  : struct with fields
%                   .coords    [2×P] int16|int32  (x;y) concatenated
%                   .offsets   [N+1×1] uint64     0-based prefix sums, N=T*S
%                   .grid_size [1×2]  uint32      [T S]
%                 (Optionally .offsets_2d, .lengths_2d; not required)
%
% Name–Value
%   'dtype'      : output dtype for values (default: 'single')
%
% Outputs
%   pixPacked : [1×P] dtype   concatenated pixel values
%               
%   averageTS : [T×S] double, mean of each segment (NaN if empty)
%
% Slicing a segment (t,s):
%   k  = (t-1)*S + s;
%   a  = double(pixPacked.offsets(k)) + 1;
%   b  = double(pixPacked.offsets(k+1));
%   v  = pixPacked.values(a:b);   % 1×K

%if nargin < 3, dtype = 'single'; end

ip = inputParser;
ip.addParameter('dtype','single',@(s)ischar(s)||isstring(s));
ip.parse(varargin{:});
outClass = char(ip.Results.dtype);

    T = double(packedModelstructure.grid_size(1));
    S = double(packedModelstructure.grid_size(2));
    N = T * S;

    % Basic sanity
    assert(numel(Chn) >= T, 'Chn must contain at least T frames.');

    frameSize = size(Chn{1});  % [H W]
P = double(packedModelstructure.offsets(end));   % total points    
pixValues = zeros(1, P, outClass);
averageTS = NaN(T, S);

    % Parallelize across sweeps (better cache reuse of Chn{t} inside inner loop)
for t = 1:T
    F = Chn{t};                     % [H×W]
    [H,W] = size(F);

    for s = 1:S
        k = (t-1)*S + s;

        a = double(packedModelstructure.offsets(k))   + 1;
        b = double(packedModelstructure.offsets(k+1));
        if b < a, 
            continue; 
        end

        % coords for this segment (2×K)
        xy = packedModelstructure.coords(:, a:b);
        x  = double(xy(1,:)).';     % K×1
        y  = double(xy(2,:)).';     % K×1

        % guard against out-of-bounds coords (should be rare after cleaning)
        in = (x>=1 & x<=W & y>=1 & y<=H);
        K  = numel(x);
        vals = NaN(K,1);

        if any(in)
            idx = sub2ind([H W], y(in), x(in));
            vals(in) = F(idx);
        end

        pixValues(1,a:b) = cast(vals, outClass);
        averageTS(t,s)   = mean(vals, 'omitnan');
    end
end
packedModelstructure.values = pixValues;
%%
if nargout > 2
    LegacyCellVals = cell(double(T), double(S));
    pp = ModelSlicingHelpers(packedModelstructure);

    for t = 1:double(T)
        for s = 1:double(S)
            v = pp.slice_values_ts(t, s);
            LegacyCellVals{t,s} = v;
            %{
            k = (t-1)*S + s;
            a = double(offsets(k))   + 1;
            b = double(offsets(k+1));
            if a <= b
                legacyCell{t,s} = [ double(coords(1,a:b)).' , double(coords(2,a:b)).' ];
            else
                legacyCell{t,s} = cast([], coordClass);
            end
            %}
        end
    end
end
end


%{
[numFrames, numSweeps] = size(outputwodups);
outputpixvals = cell(numFrames, numSweeps);
averageint = NaN(numFrames, numSweeps);
frameSize = size(Chn{1});

parfor i = 1:numSweeps
    localPixVals = cell(numFrames, 1);
    localAvg = NaN(numFrames, 1);

    for ii = 1:numFrames
        coords = outputwodups{ii, i};
        if ~isempty(coords) && size(coords, 2) == 2
            ind = sub2ind(frameSize, coords(:, 2), coords(:, 1));
            vals = Chn{ii}(ind);

            localPixVals{ii} = vals;
            localAvg(ii) = mean(vals, 'omitnan');
        else
            localPixVals{ii} = NaN;
        end
    end

    outputpixvals(:, i) = localPixVals;
    averageint(:, i) = localAvg;
end
%}
%{
profile on
pardowntime = size(outputwodups,1);
for i = 1:size(outputwodups, 2) % Sweep across vessels
    for ii = 1:pardowntime % Down time
        crosspixcoords = outputwodups{ii, i}; % Get coordinates to sample then average
        siz = size(crosspixcoords, 1);
        siz2 = size(crosspixcoords,2);
        if siz > 0 & siz2 > 1 % Only if there are coordinates
            a = sub2ind(size(Chn{1},[1,2]), crosspixcoords(:, 2), crosspixcoords(:, 1));
            outputpixvals{ii, i} = Chn{ii}(a);
            averageint(ii, i) = mean(a, 'omitnan'); % Average, omitting NaNs
        else
            % Handle the case where there are no coordinates
            outputpixvals{ii, i} = NaN; % Fill with NaNs if no coordinates
            averageint(ii, i) = NaN; % Set average to NaN
        end
    end
end
profile viewer
%}