function [BVesFilled] = FillBinaries(BVes)
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
            se = strel("line",15,0);
            se2 = strel("line",20,0);
            BVes(:,1,:) = 1;
            BVes(:,end,:) = 1;
            parfor i = 1 : size(BVes,3)
                Vesobj = imdilate(BVes(:,:,i),se);
                Vesobj = bwareafilt(Vesobj,2);
                Vesobj = imfill(Vesobj,'holes');
                Vesobj = imerode(Vesobj,se2);
                Vesobj(:,1) = 1;
                Vesobj(:,end) = 1;
                Vesobj = imfill(Vesobj,'holes');
                Vesobj(:,1) = 0;
                Vesobj(:,end) = 0;
                BVesFilled(:,:,i) = bwareafilt(Vesobj,1);
            end
            BVesFilled = smoothdata(BVesFilled,3);
            BVesFilled = imbinarize(BVesFilled,0.4);
            BVesFilled(:,1,:) = 0;
            BVesFilled(:,end,:) = 0;
end
