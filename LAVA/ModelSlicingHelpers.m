classdef ModelSlicingHelpers
    %%% USAGE:

    %pp = ModelSlicingHelpers(LAVA.Endfoot.Mask.SlicedCoordinatesAndValues)
    % pp = PackagedModel(PackagedStructure);          % wrap your existing struct
    % [a,b] = pp.seg_range(t, s);
    % v     = pp.slice_values_ts(t, s);
    % C1    = pp.slice_values_frame(t);
    % C2    = pp.slice_values_sweep(s);

    properties (SetAccess = immutable)
        coords
        grid_size   (1,2) double {mustBePositive, mustBeInteger} = [1,1];
        offsets     (:,1) uint64
        values      (1,:) % keep your chosen dtype (single/double/int)
    end

    methods
        function obj = ModelSlicingHelpers(packed)
            obj.grid_size = double(packed.grid_size);
            obj.offsets   = packed.offsets(:);     % ensure column
            obj.values    = packed.values;
            obj.coords    = packed.coords;
            T = obj.grid_size(1);
            S = obj.grid_size(2);
            assert(numel(obj.offsets) == T*S + 1, 'offsets must be T*S+1 long');
            assert(obj.offsets(1) == 0, 'offsets(1) should be 0');
            last = double(obj.offsets(end));
            assert(last <= numel(obj.values), 'offsets(end) exceeds values length');
        end
        function [a,b] = seg_range(obj, t, s)
            % 1-based [a,b] range for segment (t,s)
            T = double(obj.grid_size(1));
            S = double(obj.grid_size(2));
            k = (t-1)*S + s;
            a = double(obj.offsets(k))   + 1;
            b = double(obj.offsets(k+1));
        end

%%% RECOVER FULL MODEL TENSORS vvvvvvvvvvv
        function [V, mask] = values_tensor(obj, padval)
            % Dense values tensor shaped L×S×T (L = max segment length)
            % padval default: NaN for float types, 0 for integer types
            T = obj.grid_size(1); S = obj.grid_size(2);
            segLens = diff(double(obj.offsets));      % length T*S
            L = max(segLens);
            if nargin < 2
                if isfloat(obj.values), padval = NaN; else, padval = 0; end
            end
            V    = repmat(cast(padval, 'like', obj.values), [L, S, T]);
            mask = false(L, S, T);

            for t = 1:T
                for s = 1:S
                    [a,b] = obj.seg_range(t,s);
                    n = b - a + 1;
                    if n > 0
                        % write along the first (length) dimension
                        V(1:n, s, t)    = obj.values(1, a:b).';  % ensure column
                        mask(1:n, s, t) = true;
                    end
                end
            end
        end

        function [C4, mask] = coords_tensor(obj, padval)
            % Dense coords tensor shaped D×L×S×T (matches L×S×T values)
            % padval default: NaN for float coords, 0 otherwise
            if isempty(obj.coords)
                C4 = []; mask = []; return;
            end
            T = obj.grid_size(1); S = obj.grid_size(2);
            D = size(obj.coords,1);
            segLens = diff(double(obj.offsets)); L = max(segLens);
            if nargin < 2
                if isfloat(obj.coords), padval = NaN; else, padval = 0; end
            end
            likeC = cast(padval, 'like', obj.coords);
            C4    = repmat(likeC, [D, L, S, T]);
            mask  = false(L, S, T);

            for t = 1:T
                for s = 1:S
                    [a,b] = obj.seg_range(t,s);
                    n = b - a + 1;
                    if n > 0
                        C4(:, 1:n, s, t) = obj.coords(:, a:b);
                        mask(1:n, s, t)  = true;
                    end
                end
            end
        end

        function U = unique_coords_per_time(obj, C4)  
            if nargin < 2 || isempty(C4)
                [C4, ~] = obj.coords_tensor();  % D×L×S×T
            end
            if isempty(C4), U = {}; return; end
            assert(size(C4,1) == 2, 'First dim of coords tensor must be 2 (x;y).');

            T = size(C4,4);
            U = cell(1,T);
            for t = 1:T
                XY = reshape(C4(:,:,:,t), 2, []).';
                XY = XY(all(isfinite(XY),2), :);
                Thepixels = unique(XY, 'rows', 'stable');

                Thepixels(1,:) = [];
                %[rex,~] = find(Thepixels(:,2) > dims(1));
                %Thepixels(rex,:) = [];   %nulls out of image bounds pixels
                %[rex,~] = find(Thepixels(:,1) > dims(2));
                %Thepixels(rex,:) = [];   %nulls out of image bounds pixels
                if size(isnan(Thepixels),1) > 0
                    excludeam = isnan(Thepixels);
                    excludeam = sum(excludeam,2);
                    Thepixels(excludeam > 0,:) = [];
                end
                U{t} = Thepixels;
            end
        end

%%% SPECIFIC SLICING vvvvvvvvvvvvvv

        function v = slice_values_ts(obj, t, s)
            [a,b] = obj.seg_range(t, s);
            if b < a
                v = obj.values(1,1:0);
            else
                v = obj.values(1,a:b);
            end
        end

        function C = slice_values_frame(obj, t)
            % 1×S cell, each cell = values for (t,s)
            S = double(obj.grid_size(2));
            C = cell(1,S);
            for s = 1:S
                C{s} = obj.slice_values_ts(t, s);
            end
        end

        function C = slice_values_sweep(obj, s)
            % T×1 cell, each cell = values for (t,s)
            T = double(obj.grid_size(1));
            C = cell(T,1);
            for t = 1:T
                C{t} = obj.slice_values_ts(t, s);
            end
        end

        %Now helpers to slice coordinates into cell arrays
        function Cseg = slice_coords_ts(obj, t, s)
            % Returns D×N coords for segment (t,s)
            % If coords are empty, returns [] (or empty D×0 if defined)
            [a,b] = obj.seg_range(t, s);
            if isempty(obj.coords)
                Cseg = [];
                return;
            end
            if b < a
                % preserve dimensionality (D×0)
                Cseg = obj.coords(:, 1:0);
            else
                Cseg = obj.coords(:, a:b);
            end
        end

        function C = slice_coords_frame(obj, t)
            % 1×S cell, each cell = D×Ni coords for (t,s)
            S = double(obj.grid_size(2));
            C = cell(1,S);
            for s = 1:S
                C{s} = obj.slice_coords_ts(t, s);
            end
        end

        function C = slice_coords_sweep(obj, s)
            % T×1 cell, each cell = D×Ni coords for (t,s)
            T = double(obj.grid_size(1));
            C = cell(T,1);
            for t = 1:T
                C{t} = obj.slice_coords_ts(t, s);
            end
        end

    end
end