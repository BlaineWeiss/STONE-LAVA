function [LinearChnSample,validind, OutsideVesselScan] = SampleOverBinary(Chn1,coordy,coordx,ModelMask)
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
validind = size(Chn1,1) >= coordy' & coordy' >= 1 & size(Chn1,2) >= coordx' & coordx' >= 1;
indexes = sub2ind([size(Chn1(:,:,1))], floor(coordy'), floor(coordx'));
parfor i = 1:size(Chn1,3)
    m = Chn1(:,:,i);
   c1D = single(zeros(size(coordy')));
    c1D(validind) = m(indexes(validind));
    LinearChnSample(:,:,i) = c1D';
    ScanChn1 = double(LinearChnSample(:,:,i));
    ScanChn1(ModelMask(:,:,i)) = 0; %%%Excludes vessel from vessel area scan was 0 instead of nan
    OutsideVesselScan(:,:,i) = ScanChn1;
end
end
