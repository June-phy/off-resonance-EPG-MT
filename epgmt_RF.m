%-------------------------------------------------------------------------
% Function:    epgmt_RF
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
%   Simulates the effects of an RF pulse in the Extended Phase Graph (EPG)
%   framework. The RF pulse rotates the magnetization states (F+, F-, and Z)
%   by applying a rotation matrix that includes flip angle, off-resonance effects,
%   and phase cycling.
%
%   Inputs:
%     S      : EPG states (matrix with rows for F+, F-, and Z states)
%     theta  : Flip angle of the RF pulse (rad), equivalent to w1 * tp
%     beta   : Off-resonance angle (rad), equivalent to sqrt(w1^2 + dw^2) * tp
%     phi    : Azimuthal angle of the RF pulse (rad), defines phase cycling
%     Rrfb   : Relaxation rate during the RF pulse [1/s]
%     tp     : Duration of the RF pulse [s]
%
%   Outputs:
%     S      : Updated EPG states after applying the RF pulse
%     R      : RF rotation matrix used for the transformation
%
%   Note:
%     - The RF rotation matrix `R` is constructed using rotation operations
%       in combination with phase cycling and the azimuthal angle `phi`.
%     - Off-resonance effects are included in the construction of `R`.
%     - Relaxation during the RF pulse is modeled by an exponential decay term
%       applied to the longitudinal magnetization (Z state).
%-------------------------------------------------------------------------
function [S, R] = epgmt_RF(S, theta, beta, phi, Rrfb, tp)

R = [cos(beta),             -cos(theta)*sin(beta),                     sin(theta)*sin(beta);
     cos(theta)*sin(beta),   cos(beta)+(1-cos(beta))*(sin(theta))^2,  (1-cos(beta))*sin(theta)*cos(theta);
    -sin(theta)*sin(beta),  (1-cos(beta))*sin(theta)*cos(theta),       cos(beta)+(1-cos(beta))*(cos(theta))^2];

P = [1,  1i, 0;
     1, -1i, 0;
     0,  0, 1];

P_inv = [1,   1, 0;
        -1i, 1i, 0;
         0,   0, 2]/2;

Rp = [exp(1i*phi),  0, 0;
     0, exp(-1i*phi), 0;
     0,          0, 1];

Rp_ = [exp(-1i*phi), 0, 0;
      0,  exp(1i*phi), 0;
      0,          0, 1];

R = Rp*P*R*P_inv*Rp_;

R(4,4) = exp(-Rrfb*(tp));

S = R * S; 