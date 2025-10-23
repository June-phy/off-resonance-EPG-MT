%-------------------------------------------------------------------------
% Function:    epgmt_Ex_Relax
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
%   Simulates the effects of exchange and relaxation in the Extended Phase
%   Graph (EPG) model for two coupled pools (e.g., free and bound pools in
%   magnetization transfer).
%
%   Inputs:
%     S      : EPG state vector (4x1 matrix for F+, F-, Za, Zb)
%     R1a    : Longitudinal relaxation rate of pool a [1/s]
%     R2a    : Transverse relaxation rate of pool a [1/s]
%     R1b    : Longitudinal relaxation rate of pool b [1/s]
%     kab    : Exchange rate from pool a to pool b [1/s]
%     kba    : Exchange rate from pool b to pool a [1/s]
%     tp     : Pulse duration [s]
%     td     : Delay duration [s]
%     M0a    : Equilibrium magnetization of pool a
%     M0b    : Equilibrium magnetization of pool b
%
%   Outputs:
%     S_new  : Updated EPG state vector after relaxation and exchange
%     E      : Transition matrix for relaxation and exchange
%     B      : Offset vector due to equilibrium magnetization
%
%   Note:
%     - This function calculates the effects of relaxation and exchange
%       based on the Bloch-McConnell equations.
%     - The magnetization evolution is represented as:
%         S_new = E * S + B
%     - The input `S` should be a 4x1 vector representing the current state
%       of the EPG system.
%
%   References:
%     Please cite relevant literature on EPG modeling and Bloch-McConnell
%     equations when using this code.
%
%-------------------------------------------------------------------------
function [S_new, E, B] = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b)

T = diag([-R2a -R2a]);
Texpm = expm( (tp+td)*T ); 
L = [ [-R1a-kab kba] ; [kab -R1b-kba] ];
Lexpm = expm( (tp+td)*L);
E = blkdiag(Texpm, Lexpm);

C_L = [M0a*R1a; M0b*R1b];
Zoff = (Lexpm - eye(2))*(L\C_L);
B = zeros([4 1]);
B([3 4]) = Zoff; 

S_new = E*S + B;