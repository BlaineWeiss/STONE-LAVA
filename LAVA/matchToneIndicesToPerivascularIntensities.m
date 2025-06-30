function [normca, resizedcrosssections] = matchToneIndicesToPerivascularIntensities(Parallelregion, areasections, icarea)
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
%{

normca = normalize(Parallelregion, 2);
    desired_size = max(icarea);
    resizedcrosssections = zeros(size(Parallelregion, 2), desired_size);
    N = size(areasections, 1);
    for i = 1:N
        n = length(areasections{i});
        resizedcrosssections(i,:) = interp1(linspace(1, n, numel(areasections{i})), areasections{i}, linspace(1, n, desired_size));
    end
end
%}
%%%%%%%%%%%%%%%%%%%%%%%

profile on
    normca = normalize(Parallelregion, 2);
    desired_size = max(icarea);
    resizedcrosssections = NaN(size(Parallelregion, 2), desired_size);
    N = size(areasections, 1);
    parfor i = 1:N
        section = areasections{i};
        if isempty(section)
            continue;
        end
        resizedcrosssections(i,:) = interp1(linspace(1, numel(section), numel(section)), section, linspace(1, numel(section), desired_size), 'linear', 'extrap');
    end

profile viewer
end
