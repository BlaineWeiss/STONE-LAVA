function [Dialine2] = SampleOverModel(yi,xi,Chn2r)
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

fimax = size(xi,1)
samps = size(xi,2)
ChnLength = size(Chn2r,3);
Dialine2 = single(zeros(samps,fimax,ChnLength));
tic
Chn22 = Chn2r(:,:,1);
validind = size(Chn22,1) >= yi & yi >= 1 & size(Chn22,2) >= xi & xi >= 1;
yi(~validind) = NaN;
xi(~validind) = NaN;
try
    indexes = sub2ind([size(Chn2r,[1,2])], yi, xi); %LINEARLY INDEXES LINE ON WHOLE IMAGE  (add condition to shrink x if subscript out of range)
catch
    for fi = 1:fimax
        try
            indexes = sub2ind([size(Chn2r,[1,2])], yi(fi,:), xi(fi,:));
        end
    end
end
parfor stack = 1: ChnLength
    Dialine = single(zeros(fimax,samps));
    Chn22 = Chn2r(:,:,stack);
    %m = zeros(512,512);
    m = Chn22; %(floor(NDiam.YData),floor(NDiam.XData)) %(subimage of line field)
    Dialine(validind) = m(indexes(validind));
    Dialine2(:,:,stack) = Dialine';
end
end
