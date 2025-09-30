classdef ModelSlicingHelpers
%%% USAGE:

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

        function v = slice_values_ts(obj, t, s)
            [a,b] = obj.seg_range(t, s);
            if b < a 
                v = obj.values(1,1:0); 
            else 
                v = obj.values(1,a:b); 
            end
        end

        function C = slice_values_frame(pixPacked, t)
            % 1×S cell, each cell = values for (t,s)
            S = double(obj.grid_size(2));
            C = cell(1,S);
            for s = 1:S 
                C{s} = obj.slice_values_ts(t, s);
            end
        end

        function C = slice_values_sweep(pixPacked, s)
            % T×1 cell, each cell = values for (t,s)
            T = double(obj.grid_size(1));
            C = cell(T,1);
            for t = 1:T 
                C{t} = obj.slice_values_ts(t, s); 
            end
        end

    end
end