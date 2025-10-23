%-------------------------------------------------------------------------
% Function:    epgmt_Grad
% Authors:     Qianxue Shan
% Email:       qxshan@link.cuhk.edu.hk
% Version:     1.0
% Date:        2025-10-23
%
% Copyright (c) 2025 Qianxue Shan. All rights reserved.
%
% License and Usage Notice:
%   This code is provided strictly for academic and research purposes only.
%   Any commercial use, including but not limited to sale, redistribution,
%   or integration into proprietary software, is strictly prohibited without
%   explicit written permission from the authors.
%
%   Modification of this code, its header comments, or removal of this notice,
%   in whole or in part, is EXPRESSLY FORBIDDEN without prior written consent
%   from the authors.
%
%   By using, copying, or referencing this code, you agree to abide by these terms.
%   For any inquiries or requests, please contact the authors at the email above.
%
% Description:
%   Simulates the effects of a gradient pulse in the Extended Phase Graph (EPG)
%   framework by updating the EPG states. The gradient shifts the phase states
%   (F+ and F-) and ensures proper handling of boundary conditions.
%
%   Inputs:
%     S      : EPG states (matrix with rows for F+, F-, Za, and Zb)
%
%   Outputs:
%     S      : Updated EPG states after applying the gradient pulse
%
%   Note:
%     - The gradient pulse updates the F+ and F- states by shifting them to
%       higher and lower phase states, respectively.
%     - Boundary conditions are handled such that the lowest F+ state is
%       conjugated with the corresponding F- state.
%
%-------------------------------------------------------------------------
function S = epgmt_Grad(S)

S = [S [0;0;0;0]];
S(1,:) = circshift(S(1,:),[0 1]);	% Shift Fp states.
S(2,:) = circshift(S(2,:),[0 -1]);	% Shift Fm states.
S(2,end)=0;					        % Zero highest Fm state.
S(1,1) = conj(S(2,1));			    % Fill in lowest Fp state.