function [sumend, Basefluoroend, ChangeFluoroend, ratiofluoroend] = calculateFluorescenceChanges(EndfeetonChannelMask, stimphase)
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
Endfeet2 = EndfeetonChannelMask;
Endfeet2(Endfeet2 == 0) = NaN;
sumend = mean(Endfeet2, 1, 'omitnan');
Basefluoroend = mean(sumend(1:stimphase), 'omitnan');
ChangeFluoroend = (sumend - Basefluoroend) ./ Basefluoroend;
ratiofluoroend = 100 .* (sumend ./ Basefluoroend);
end
