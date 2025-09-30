function [packed, legacyCell] = extractUniqueROICoords(ROIcoordx, ROIcoordy, varargin)

%function outputwodups = extractUniqueROICoords(ROIcoordx, ROIcoordy)
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

% Refactor: produce a compact "packed" representation of SlicedCoordinates
% that is fast to write/read and HDF5-friendly.
%
% INPUT
%   ROIcoordx, ROIcoordy : [numPoints x numSweeps x numFrames] double
%       Points are samples along each cross-section sweep over time.
%
% OUTPUT (default)
%   packed : struct with fields
%       .coords       : [2 x P] int16|int32, concatenated [x;y] for all cells
%       .offsets      : [N+1 x 1] uint64, prefix sum; segment k is columns (offsets(k)+1 : offsets(k+1))
%       .t_index      : [N x 1] uint32, time index (row) for each segment (1..T)
%       .s_index      : [N x 1] uint32, sweep index (col) for each segment (1..S)
%       .grid_size    : [1 x 2] uint32, [T S]
%       .coords_dtype : string, 'int16' or 'int32'
%       .order        : string, 'row-major' (k = (t-1)*S + s)
%       .note         : string, brief doc
%
% OPTIONAL 2nd OUTPUT
%   legacyCell : {T x S} cell, each cell = [K x 2] int16|int32 (for compatibility)
%
% USAGE
%   packed = extractUniqueROICoords(ROIcoordx, ROIcoordy);
%   [packed, legacy] = extractUniqueROICoords(ROIcoordx, ROIcoordy);  % also returns legacy cell grid
%
% Motivation: A T×S cell of tiny [K×2] arrays causes massive per-cell overhead.
% Packing into one big numeric + offsets reduces memory & speeds up HDF5.

% ---- Parse & check dims ---------------------------------------------------
[numPoints, numSweeps, numFrames] = size(ROIcoordx);
assert(isequal(size(ROIcoordx), size(ROIcoordy)), 'ROIcoordx/ROIcoordy must match size.');

T = uint32(numFrames);  % "time" axis = rows in your old grid
S = uint32(numSweeps);  % "sweep/cross-section" axis = columns in your old grid
N = double(T) * double(S);  % total segments (cells) in row-major order

% We need the max coordinate magnitude to choose int16 vs int32
% Do a quick streamed pass to compute per-segment counts and the max |coord|
segCounts = zeros(N,1, 'uint32');   % points per segment
maxAbsCoord = 0;                    % track max abs(x/y) across all segments

% Row-major index helper (no storage): k = (t-1)*S + s, with t∈[1..T], s∈[1..S]
k = @(t,s) (double(t)-1)*double(S) + double(s);

% ------------- First pass: count points & discover max magnitude -----------
for t = 1:double(T)
    % Slice once along time for efficiency
    X_t = squeeze(ROIcoordx(:, :, t));   % [numPoints x S]
    Y_t = squeeze(ROIcoordy(:, :, t));   % [numPoints x S]

    for s = 1:double(S)
        cols = [X_t(:,s), Y_t(:,s)];

        % Original cleaning: floor, remove any row with zeros or NaNs
        cols = floor(cols);
        cols = cols(all(cols,2) & ~any(isnan(cols),2), :);

        % Unique rows, stable (preserve original order)
        if ~isempty(cols)
            cols = unique(cols, 'rows', 'stable');
            maxAbsCoord = max(maxAbsCoord, max(abs(cols),[],'all'));
            segCounts(k(t,s)) = uint32(size(cols,1));
        else
            segCounts(k(t,s)) = uint32(0);
        end
    end
end

% Decide dtype for coordinates
if maxAbsCoord <= intmax('int16')
    coordClass = 'int16';
else
    coordClass = 'int32';
end

% Total points P (sum of K for all segments)
P = double(sum(segCounts, 'native'));

% Preallocate packed arrays
coords  = zeros(2, P, coordClass);     % concatenated [x; y]
offsets = zeros(N+1, 1, 'uint64');     % prefix sums (start at 0)

% ------------- Second pass: fill coords + offsets + indices ---------------
writeHead = uint64(0);   % number of points already written into coords

for t = 1:double(T)
    X_t = squeeze(ROIcoordx(:, :, t));
    Y_t = squeeze(ROIcoordy(:, :, t));

    for s = 1:double(S)
        kk = k(t,s);            % linear segment id in row-major
        nThis = uint64(segCounts(kk));

        % Set prefix offset for this segment (offsets(kk) stores the start)
        offsets(kk) = writeHead;

        if nThis == 0
            continue
        end

        cols = [X_t(:,s), Y_t(:,s)];
        cols = floor(cols);
        cols = cols(all(cols,2) & ~any(isnan(cols),2), :);
        cols = unique(cols, 'rows', 'stable');

        % Write slice into coords: columns (writeHead+1) .. (writeHead+nThis)
        colRange = (writeHead+1) : (writeHead+nThis);
        coords(1, colRange) = cast(cols(:,1), coordClass);
        coords(2, colRange) = cast(cols(:,2), coordClass);

        writeHead = writeHead + nThis;   % advance head
    end
end

% Final sentinel offset
offsets(N+1) = writeHead;

% Sanity: writeHead should equal P
% (If not, something changed between the two passes.)
% assert(writeHead == uint64(P), 'Packed count mismatch.');

% ------------- Build packed struct ----------------------------------------
packed = struct( ...
    'coords',        coords, ...
    'offsets',       offsets, ...          % length N+1, 0-based prefix sums
    'grid_size',     uint32([T S]), ...
    'coords_dtype',  string(coordClass), ...
    'order',         "row-major", ...
    'note',          "Segment k = (t-1)*S + s; segment coords are coords(:, offsets(k)+1 : offsets(k+1))." ...
);

% ---- Optional convenience: 2-D starts/lengths (derived from 1-D offsets) ----
N = T * S;

starts0 = packed.offsets(1:end-1);   % uint64, 0-based starts, length N
lens    = diff(packed.offsets);      % uint64 lengths, length N

packed.offsets_2d = reshape(starts0, [T S]);      % uint64 [T x S], 0-based
packed.lengths_2d = reshape(uint32(lens), [T S]); % uint32 [T x S]

% (Optional) If you like 1-based starts for MATLAB slicing:
% packed.starts1_2d = packed.offsets_2d + 1;  % uint64 [T x S]





% ------------- Optional legacy cell (only if requested) ---REPLACE WITH STRUCTURE CLASS METHOD CALLS..
if nargout > 1
    legacyCell = cell(double(T), double(S));
    for t = 1:double(T)
        for s = 1:double(S)
            kk = k(t,s);
            a = double(offsets(kk))   + 1;
            b = double(offsets(kk+1));
            if a <= b
                legacyCell{t,s} = [ double(coords(1,a:b)).' , double(coords(2,a:b)).' ];
            else
                legacyCell{t,s} = cast([], coordClass);
            end
        end
    end
end
end








%%
%{
[numPoints, numSweeps, numFrames] = size(ROIcoordx);
outputwodups = cell(numFrames, numSweeps);
for i = 1:numSweeps
    Sweepix = squeeze(ROIcoordx(:, i, :));
    Sweepiy = squeeze(ROIcoordy(:, i, :));

    for sta = 1:numFrames
        aa = floor([Sweepix(:, sta), Sweepiy(:, sta)]);
        aa = aa(all(aa, 2) & ~any(isnan(aa), 2), :);


        % coordinate list becomes cell arrays inside 2d array struct.. Not good for hdf5 creation
        % Get 
        if ~isempty(aa)
            outputwodups{sta, i} = int16(unique(aa, 'rows', 'stable'));
        else
            outputwodups{sta, i} = int16([]);
        end
    end
end
%}