%-------------------------------------------------------------------------
% Script:      EPG DANTE MT Simulation - Liver vs Blood
% Authors:     Qianxue Shan, Weitian Chen
% Email:       qxshan@link.cuhk.edu.hk
% Version:     1.0
% Date:        2025-10-23
%
% Copyright (c) 2025 Qianxue Shan. All rights reserved.
%
% License and Usage Notice:
%   This script is provided strictly for academic and research purposes only.
%   Any commercial use, including but not limited to sale, redistribution,
%   or integration into proprietary software, is strictly prohibited without
%   explicit written permission from the authors.
%
%   Modification of this script, its header comments, or removal of this notice,
%   in whole or in part, is EXPRESSLY FORBIDDEN without prior written consent
%   from the authors.
%
%   By using, copying, or referencing this script, you agree to abide by these terms.
%   For any inquiries or requests, please contact the authors at the email above.
%
% Description:
%   This script simulates magnetization transfer (MT) effects for liver and blood 
%   tissues using DANTE pulse trains, based on the Extended Phase Graph (EPG) 
%   framework and Bloch-McConnell equations.
%
%   Main Features:
%     - Simulates the effects of varying flow angles, gradient amplitudes, 
%       pulse numbers, and B1 inhomogeneity on tissue magnetization.
%     - Produces a 2×4 figure comparing magnetization evolution for liver (row 1)
%       and blood (row 2) under different experimental conditions.
%     - Includes local functions for running angle sweep, gradient sweep, 
%       pulse number evolution, and B1 sweep simulations.
%
%-------------------------------------------------------------------------
%%
clear; close all; clc;

%% ---------------------- Common sequence settings ------------------------
gamma      = 2.6752e8;             % rad/(s*T)
tp         = 0.7e-3;               % RF pulse duration [s]
td         = 1.0e-3;               % spacing between DANTE modules [s]
tg         = 1e-3;                 % gradient duration [s]
phi        = 0;                    % RF phase [rad]
dw         = 2*pi*1000;            % off-resonance [rad/s]
w1_base    = 2*pi*100;             % base RF nutation rate [rad/s]
Gz_nominal = 5e-3;                 % nominal gradient amplitude [T/m]

% Flow velocities (cm/s for legend; converted to m/s for computation)
v_arr_cm = [0, 0.1, 0.5, 1.0, 5.0, 10.0];
v_arr    = v_arr_cm * 1e-2;        % m/s
n_v      = numel(v_arr);

% Angle sweep (deg & rad)
angle_deg = linspace(0, 180, 181);
angle_rad = deg2rad(angle_deg);

% Gradient sweep (mT/m for display, T/m for computation)
Gz_mTm = linspace(0, 50, 101);
Gz_arr = Gz_mTm * 1e-3;

% B1 scaling factors
b1_scale = 0.6:0.1:1.4;
Np_b1    = 200;

% Np evolution (number of pulses)
Np_max  = 1000;
Np_axis = 0:Np_max;

%% ---------------------- Tissue parameter definitions --------------------
tissues = {
    struct('name', 'Liver', ...
           'T1',   0.812, ...            % [s]
           'T2',   0.042, ...            % [s]
           'T2b',  7.7e-6, ...           % [s]
           'fb',   0.069, ...
           'kba',  51)
    struct('name', 'Blood', ...
           'T1',   1932*1e-3, ...        % [s]
           'T2',   275*1e-3, ...         % [s]
           'T2b',  280*1e-6, ...         % [s]
           'fb',   0.028, ...
           'kba',  35)
};

%% ---------------------- Prepare figure & layout -------------------------
bigFig = figure('Color', 'w', ...
                'Name', 'EPG DANTE MT - Liver (row 1) vs Blood (row 2)', ...
                'Position', [100, 100, 1400, 700]);
tl = tiledlayout(bigFig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Labels for subplots (a)-(h)
tile_labels = {'(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)'};

%% ---------------------- Loop over tissues & generate subplots -----------
for idx = 1:numel(tissues)
    tissue = tissues{idx};
    fprintf('Processing tissue: %s\n', tissue.name);

    % Derived MT relaxation parameters
    T1  = tissue.T1;
    T2  = tissue.T2;
    T2b = tissue.T2b;
    fb  = tissue.fb;
    kba = tissue.kba;
    kab = fb * kba;

    R1a = 1 / T1;
    R2a = 1 / T2;
    R1b = R1a;

    M0a = 1;
    M0b = fb;

    % Precompute RF parameters for baseline (used in angle/gradient/Np scans)
    [theta_base, beta_base, Rrfb_base] = prep_rf_params(w1_base, dw, T2b, tp);

    % Common settings
    Np_angle      = 200;
    Np_grad       = 200;
    angle_fix_rad = 0;  % flow parallel to +z gradient

    % Determine tile offset (row-major ordering)
    tile_offset = (idx - 1) * 4;

    %% 1) Angle sweep subplot (column 1)
    ax = nexttile(tl, tile_offset + 1);
    hold(ax, 'on'); grid(ax, 'on');

    for j = 1:n_v
        v = v_arr(j);
        Mz_angle = run_angle_sweep(Np_angle, angle_rad, v, ...
                                   tp, td, tg, ...
                                   R1a, R2a, R1b, kab, kba, ...
                                   M0a, M0b, ...
                                   theta_base, beta_base, phi, Rrfb_base, ...
                                   Gz_nominal);
        plot(ax, angle_deg, Mz_angle, 'LineWidth', 1.4);
    end

    if idx == 2
        xlabel(ax, '\theta (degrees)');
        ax.XLabel.FontSize = 14;
    else
        xlabel(ax, '');
    end

    ylabel(ax, 'M_z');
    ax.YLabel.FontSize = 14;

    title(ax, sprintf('%s %s - Flow Direction', ...
                      tile_labels{tile_offset + 1}, tissue.name));
    ax.Title.FontSize   = 14;
    ax.Title.FontWeight = 'bold';
    if idx == 1
        leg = legend(ax, arrayfun(@(vc) sprintf('v = %.1f cm/s', vc), v_arr_cm, ...
                                  'UniformOutput', false), 'Location', 'southwest', 'Box', 'off');
        leg.FontSize = 12;
    end

    ax.XLim  = [0, 180];
    ax.XTick = 0:30:180;
    ax.FontSize = 12;

    %% 2) Gradient sweep subplot (column 2)
    ax = nexttile(tl, tile_offset + 2);
    hold(ax, 'on'); grid(ax, 'on');

    for j = 1:n_v
        v = v_arr(j);
        Mz_grad = run_gradient_sweep(Np_grad, Gz_arr, angle_fix_rad, v, ...
                                     tp, td, tg, ...
                                     R1a, R2a, R1b, kab, kba, ...
                                     M0a, M0b, ...
                                     theta_base, beta_base, phi, Rrfb_base);
        plot(ax, Gz_mTm, Mz_grad, 'LineWidth', 1.4);
    end

    if idx == 2
        xlabel(ax, 'G_z (mT/m)');
        ax.XLabel.FontSize = 14;
    else
        xlabel(ax, '');
    end

    ylabel(ax, '');
    title(ax, sprintf('%s %s - Gradient', ...
                      tile_labels{tile_offset + 2}, tissue.name));
    ax.Title.FontSize   = 14;
    ax.Title.FontWeight = 'bold';

    ax.XLim = [min(Gz_mTm), max(Gz_mTm)];
    ax.FontSize = 12;

    %% 3) Np evolution subplot (column 3)
    ax = nexttile(tl, tile_offset + 3);
    hold(ax, 'on'); grid(ax, 'on');

    for j = 1:n_v
        v = v_arr(j);
        Mz_Np = run_Np_evolution(Np_max, angle_fix_rad, Gz_nominal, v, ...
                                 tp, td, tg, ...
                                 R1a, R2a, R1b, kab, kba, ...
                                 M0a, M0b, ...
                                 theta_base, beta_base, phi, Rrfb_base);
        plot(ax, Np_axis, Mz_Np, 'LineWidth', 1.2);
    end

    if idx == 2
        xlabel(ax, 'N_p');
        ax.XLabel.FontSize = 14;
    else
        xlabel(ax, '');
    end

    ylabel(ax, '');
    title(ax, sprintf('%s %s - Number of Pulse', ...
                      tile_labels{tile_offset + 3}, tissue.name));
    ax.Title.FontSize   = 14;
    ax.Title.FontWeight = 'bold';

    ax.XLim = [0, Np_max];
    ax.FontSize = 12;

    %% 4) B1 sweep subplot (column 4)
    ax = nexttile(tl, tile_offset + 4);
    hold(ax, 'on'); grid(ax, 'on');

    for j = 1:n_v
        v = v_arr(j);
        [Mz_b1, B1_uT] = run_b1_sweep(Np_b1, b1_scale, angle_fix_rad, Gz_nominal, v, ...
                                      tp, td, tg, gamma, ...
                                      R1a, R2a, R1b, kab, kba, ...
                                      M0a, M0b, ...
                                      w1_base, dw, phi, T2b);
        plot(ax, b1_scale, Mz_b1, 'LineWidth', 1.4);
    end

    if idx == 2
        xlabel(ax, 'B_1 scale');
        ax.XLabel.FontSize = 14;
    else
        xlabel(ax, '');
    end

    ylabel(ax, '');
    title(ax, sprintf('%s %s - B_1 inhomogeneity', ...
                      tile_labels{tile_offset + 4}, tissue.name));
    ax.Title.FontSize   = 14;
    ax.Title.FontWeight = 'bold';

    ax.XLim = [min(b1_scale), max(b1_scale)];
    ax.FontSize = 12;
end

fprintf('All subplots generated in a single 2×4 figure.\n');

%% ========================================================================
%                               Local functions
% ========================================================================

function Mz_vs_angle = run_angle_sweep(Np, angle_rad, v, ...
                                       tp, td, tg, ...
                                       R1a, R2a, R1b, kab, kba, ...
                                       M0a, M0b, ...
                                       theta, beta, phi, Rrfb, ...
                                       Gz_nominal)
    n_angle = numel(angle_rad);
    Mz_vs_angle = zeros(1, n_angle);

    for idx = 1:n_angle
        angle = angle_rad(idx);
        S = [0; 0; M0a; M0b];

        for n = 1:Np
            S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
            S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
            S = epgmt_Grad(S);
            S = epgmt_Flow(S, Gz_nominal, tg, v, angle);
        end

        Mz_vs_angle(idx) = real(S(3));
    end
end

function Mz_vs_Gz = run_gradient_sweep(Np, Gz_arr, angle_fix_rad, v, ...
                                       tp, td, tg, ...
                                       R1a, R2a, R1b, kab, kba, ...
                                       M0a, M0b, ...
                                       theta, beta, phi, Rrfb)
    n_Gz = numel(Gz_arr);
    Mz_vs_Gz = zeros(1, n_Gz);

    for idx = 1:n_Gz
        Gz = Gz_arr(idx);
        S = [0; 0; M0a; M0b];

        for n = 1:Np
            S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
            S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
            S = epgmt_Grad(S);
            S = epgmt_Flow(S, Gz, tg, v, angle_fix_rad);
        end

        Mz_vs_Gz(idx) = real(S(3));
    end
end

function Mz_evolution = run_Np_evolution(Np_max, angle_fix_rad, Gz_nominal, v, ...
                                         tp, td, tg, ...
                                         R1a, R2a, R1b, kab, kba, ...
                                         M0a, M0b, ...
                                         theta, beta, phi, Rrfb)
    Mz_evolution = zeros(1, Np_max + 1);
    S = [0; 0; M0a; M0b];
    Mz_evolution(1) = real(S(3));

    for n = 1:Np_max
        S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
        S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
        S = epgmt_Grad(S);
        S = epgmt_Flow(S, Gz_nominal, tg, v, angle_fix_rad);

        Mz_evolution(n + 1) = real(S(3));
    end
end

function [Mz_vs_b1, B1_uT] = run_b1_sweep(Np, b1_scale, angle_fix_rad, Gz_nominal, v, ...
                                          tp, td, tg, gamma, ...
                                          R1a, R2a, R1b, kab, kba, ...
                                          M0a, M0b, ...
                                          w1_base, dw, phi, T2b)
    n_b1 = numel(b1_scale);
    Mz_vs_b1 = zeros(1, n_b1);
    B1_uT    = (w1_base * b1_scale / gamma) * 1e6;

    for idx = 1:n_b1
        w1 = w1_base * b1_scale(idx);
        [theta, beta, Rrfb] = prep_rf_params(w1, dw, T2b, tp);

        S = [0; 0; M0a; M0b];
        for n = 1:Np
            S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
            S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
            S = epgmt_Grad(S);
            S = epgmt_Flow(S, Gz_nominal, tg, v, angle_fix_rad);
        end

        Mz_vs_b1(idx) = real(S(3));
    end
end

function [theta, beta, Rrfb] = prep_rf_params(w1, dw, T2b, tp)
    Omega_eff = sqrt(w1^2 + dw^2);
    beta      = Omega_eff * tp;
    cos_theta = max(min(dw / Omega_eff, 1), -1);  % numerical safety
    theta     = acos(cos_theta);
    Rrfb      = RF_MT(T2b, w1, dw, 'SuperLorentzian');
end