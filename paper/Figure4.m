%-------------------------------------------------------------------------
% Script:      Demo script for Off-resonance EPG-MT Algorithm
% Authors:     Qianxue Shan, Weitian Chen
% Email:       qxshan@link.cuhk.edu.hk
% Date:        2026-01-27
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
%-------------------------------------------------------------------------
%%
clear; close all; clc;

%% ==================== USER SETTINGS ====================
% Toggle for bold labels and captions
% Set to true for bold text, false for normal text
USE_BOLD_TEXT = true;  % <-- Change this to false if you don't want bold text

% Line colors for different velocities
line_colors = [0.0000, 0.4470, 0.7410;   % 蓝色
               0.8500, 0.3250, 0.0980;   % 橙色
               0.9290, 0.6940, 0.1250;   % 黄色
               0.4940, 0.1840, 0.5560;   % 紫色
               0.4660, 0.6740, 0.1880;   % 绿色
               0.3010, 0.7450, 0.9330];  % 浅蓝色

% Line width
line_width = 1.4;

% Font sizes
title_font_size = 14;
label_font_size = 14;
axis_font_size  = 12;
legend_font_size = 12;
%% =======================================================

%% Set font weight based on toggle
if USE_BOLD_TEXT
    fontWeight = 'bold';
else
    fontWeight = 'normal';
end

%% ---------------------- Common sequence settings ------------------------
gamma      = 2.6752e8;             % rad/(s*T)
tp         = 0.7e-3;               % RF pulse duration [s]
td         = 2e-3;                 % spacing between DANTE modules [s]
tg         = 1e-3;                 % gradient duration [s]
phi        = 0;                    % RF phase [rad]
dw         = 2*pi*1000;            % off-resonance [rad/s]
w1_base    = 2*pi*100;             % base RF nutation rate [rad/s]
Gz_nominal = 25e-3;                % nominal gradient amplitude [T/m]

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
Np_b1    = 100;

% Np evolution (number of pulses)
Np_max  = 1000;
Np_axis = 0:Np_max;

%% ---------------------- Tissue parameter definitions --------------------
tissues = {
    struct('name', 'WM', ...
           'T1',   1084e-3, ...          % [s]
           'T2',   69e-3, ...            % [s]
           'T2b',  10e-6, ...            % [s]
           'fb',   0.139, ...
           'kba',  23)
    struct('name', 'Blood', ...
           'T1',   1932e-3, ...          % [s]
           'T2',   275e-3, ...           % [s]
           'T2b',  280e-6, ...           % [s]
           'fb',   0.028, ...
           'kba',  35)
};

%% ---------------------- Prepare figure & layout -------------------------
bigFig = figure('Color', 'w', ...
                'Name', 'EPG DANTE MT - White Matter (row 1) vs Blood (row 2)', ...
                'Position', [100, 100, 1400, 700]);
tl = tiledlayout(bigFig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Labels for subplots (a)-(h)
tile_labels = {'(A)', '(B)', '(C)', '(D)', '(E)', '(F)', '(G)', '(H)'};

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
    Np_angle      = 100;
    Np_grad       = 100;
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
        if iscell(line_colors)
            plot(ax, angle_deg, Mz_angle, 'LineWidth', line_width, 'Color', line_colors{j});
        else
            plot(ax, angle_deg, Mz_angle, 'LineWidth', line_width, 'Color', line_colors(j, :));
        end
    end

    if idx == 2
        xlabel(ax, '\theta (degrees)', 'FontWeight', fontWeight);
        ax.XLabel.FontSize = label_font_size;
    else
        xlabel(ax, '');
    end

    ylabel(ax, 'M_z', 'FontWeight', fontWeight);
    ax.YLabel.FontSize = label_font_size;

    title(ax, sprintf('%s %s - Flow Direction', ...
                      tile_labels{tile_offset + 1}, tissue.name), ...
          'FontWeight', fontWeight);
    ax.Title.FontSize = title_font_size;

    ax.XLim  = [0, 180];
    ax.XTick = 0:30:180;
    ax.FontSize = axis_font_size;

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
        if iscell(line_colors)
            plot(ax, Gz_mTm, Mz_grad, 'LineWidth', line_width, 'Color', line_colors{j});
        else
            plot(ax, Gz_mTm, Mz_grad, 'LineWidth', line_width, 'Color', line_colors(j, :));
        end
    end

    if idx == 2
        xlabel(ax, 'G_z (mT/m)', 'FontWeight', fontWeight);
        ax.XLabel.FontSize = label_font_size;
    else
        xlabel(ax, '');
    end

    ylabel(ax, '');
    title(ax, sprintf('%s %s - Gradient', ...
                      tile_labels{tile_offset + 2}, tissue.name), ...
          'FontWeight', fontWeight);
    ax.Title.FontSize = title_font_size;

    ax.XLim = [min(Gz_mTm), max(Gz_mTm)];
    ax.FontSize = axis_font_size;

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
        if iscell(line_colors)
            plot(ax, Np_axis, Mz_Np, 'LineWidth', line_width, 'Color', line_colors{j});
        else
            plot(ax, Np_axis, Mz_Np, 'LineWidth', line_width, 'Color', line_colors(j, :));
        end
    end

    if idx == 2
        xlabel(ax, 'N_p', 'FontWeight', fontWeight);
        ax.XLabel.FontSize = label_font_size;
    else
        xlabel(ax, '');
    end

    ylabel(ax, '');
    title(ax, sprintf('%s %s - Number of Pulse', ...
                      tile_labels{tile_offset + 3}, tissue.name), ...
          'FontWeight', fontWeight);
    ax.Title.FontSize = title_font_size;

    ax.XLim = [0, Np_max];
    ax.FontSize = axis_font_size;

    % Add legend to the third subplot (Np evolution) of the first row only
    if idx == 1
        leg = legend(ax, arrayfun(@(vc) sprintf('v = %.1f cm/s', vc), v_arr_cm, ...
                                  'UniformOutput', false), 'Location', 'east', 'Box', 'off');
        leg.FontSize = legend_font_size;
    end

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
        if iscell(line_colors)
            plot(ax, b1_scale, Mz_b1, 'LineWidth', line_width, 'Color', line_colors{j});
        else
            plot(ax, b1_scale, Mz_b1, 'LineWidth', line_width, 'Color', line_colors(j, :));
        end
    end

    if idx == 2
        xlabel(ax, 'B_1 scale', 'FontWeight', fontWeight);
        ax.XLabel.FontSize = label_font_size;
    else
        xlabel(ax, '');
    end

    ylabel(ax, '');
    title(ax, sprintf('%s %s - B_1 inhomogeneity', ...
                      tile_labels{tile_offset + 4}, tissue.name), ...
          'FontWeight', fontWeight);
    ax.Title.FontSize = title_font_size;

    ax.XLim = [min(b1_scale), max(b1_scale)];
    ax.FontSize = axis_font_size;
end

fprintf('All subplots generated in a single 2×4 figure.\n');

% Export figure
exportgraphics(gcf, 'figure parameter selection1.png', 'Resolution', 600);

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