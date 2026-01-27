%-------------------------------------------------------------------------
% Function:    epgmt_Flow
% Authors:     Qianxue Shan, Weitian Chen
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
%   Simulates the effects of flow in an Extended Phase Graph (EPG) framework.
%   This function applies flow-induced phase shifts to EPG states based on
%   the gradient strength, flow velocity, and angle between the flow direction
%   and the gradient.
%
%   Inputs:
%     S      : EPG states (matrix with rows for F+, F-, Za, and Zb)
%     Gz     : Gradient strength [T/m]
%     tg     : Gradient duration [s]
%     v      : Flow velocity [m/s]
%     angle  : Angle between flow direction and gradient [rad]
%
%   Outputs:
%     S      : Updated EPG states after applying flow-induced phase shifts
%
%   Note:
%     This implementation is based on the formalism described in the Extended
%     Phase Graph (EPG) model, accounting for flow effects according to 
%     phase modulation equations (e.g., Eq. 53 in relevant literature).
%
%   References:
%     Please cite the appropriate EPG flow modeling literature when using
%     this code.
%
%-------------------------------------------------------------------------
function S = epgmt_Flow(S, Gz, tg, v, angle)

gamma = 2.6752*1e8;
Ag = Gz*tg; 

k = gamma * Ag;

Findex = 0:length(S(1,:))-1;

exp_L = Findex * k * v * tg * cos(angle);      % J^L = exp(-1i * exp_L) (eq 53b)
exp_Tp = (Findex + 0.5) * k * v * tg * cos(angle); % J^T = exp(-1i * exp_T) (53a)
exp_Tm = (Findex - 0.5) * k * v * tg * cos(angle); % opposite index-value

S(1,:) = S(1,:) .* exp(-1i * exp_Tp);       % flow effects of F+ states
S(2,:) = S(2,:) .* exp(-1i * exp_Tm);       % flow effects of F- states
S(3,:) = S(3,:) .* exp(-1i * exp_L);        % flow effects of Za states
S(4,:) = S(4,:) .* exp(-1i * exp_L);        % flow effects of Zb states