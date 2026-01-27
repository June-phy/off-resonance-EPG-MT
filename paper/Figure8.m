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

clear; close all; clc;

%% ---------------------- Fixed experimental parameters -------------------
tp     = 0.7e-3;           % RF pulse width [s]
td     = 2.0e-3;           % Module spacing [s]
Gz     = 32e-3;            % Gradient amplitude [T/m]
tg     = 1e-3;             % Gradient duration [s]
phi    = 0;                % RF phase [rad]
angle  = 0;                % Flow direction (0 -> along +z)
gamma  = 2*pi*42.577e6;    % Proton gyromagnetic ratio [rad/(s*T)]

%% ---------------------- Tissue parameter sets ---------------------------
% Gray-matter–like
T1a_gm = 1820e-3;          % [s]
T2a_gm = 99e-3;            % [s]
T2b_gm = 9.1e-6;           % [s]
fb_gm  = 0.050;
kba_gm = 40;               % [s^-1]

% White-matter–like
T1a_wm = 1084e-3;          % [s]
T2a_wm = 69e-3;            % [s]
T2b_wm = 10e-6;            % [s]
fb_wm  = 0.139;
kba_wm = 23;               % [s^-1]

M0a_fixed = 1;

%% ---------------------- RF / offset settings ----------------------------
% High-power frequency and offset (500 Hz, 5000 Hz)
hi_freq_Hz_nom   = 500;    % High power nominal frequency [Hz]
hi_offset_Hz_nom = 5000;   % High power nominal offset [Hz]

% Low-power frequency and offset (100 Hz, 1000 Hz)
lo_freq_Hz_nom   = 100;    % Low power nominal frequency [Hz]
lo_offset_Hz_nom = 1000;   % Low power nominal offset [Hz]

% Convert to rad/s
w1_hi_nom = 2*pi*hi_freq_Hz_nom;
dw_hi_nom = 2*pi*hi_offset_Hz_nom;

w1_lo_nom = 2*pi*lo_freq_Hz_nom;
dw_lo_nom = 2*pi*lo_offset_Hz_nom;

%% ---------------------- B1 scaling (fixed) ------------------------------
b1_scale_fixed = 1.00;     % B1 inhomogeneity scaling
w1_hi_scaled   = w1_hi_nom * b1_scale_fixed;
dw_hi_scaled   = dw_hi_nom;
w1_lo_scaled   = w1_lo_nom * b1_scale_fixed;
dw_lo_scaled   = dw_lo_nom;

%% ---------------------- Initial magnetization signs ---------------------
initSigns = [+1, -1, +1, -1];

%% ---------------------- Fixed Np and velocity ---------------------------
Np = 100;                  % Fixed number of DANTE pulses
v  = 0;                    % Static condition (v = 0 m/s)

%% ---------------------- Parameter grids ---------------------------------
T1a_values = linspace(500e-3, 2500e-3, 100);   % T1a grid [s]
T2a_values = linspace(50e-3, 120e-3, 100);      % T2a grid [s]
kba_values = linspace(0, 60, 100);           % kba grid [s^-1]
T2b_values = linspace(5e-6, 15e-6, 100);       % T2b grid [s]
fb_values  = linspace(0, 0.15, 100);        % macromolecular pool fraction grid

%% ---------------------- Result containers --------------------------------
Rmpfsl_T1a_gm = NaN(size(T1a_values));
Rmpfsl_T2a_gm = NaN(size(T2a_values));
Rmpfsl_kba_gm = NaN(size(kba_values));
Rmpfsl_T2b_gm = NaN(size(T2b_values));
Rmpfsl_fb_gm  = NaN(size(fb_values));

Rmpfsl_T1a_wm = NaN(size(T1a_values));
Rmpfsl_T2a_wm = NaN(size(T2a_values));
Rmpfsl_kba_wm = NaN(size(kba_values));
Rmpfsl_T2b_wm = NaN(size(T2b_values));
Rmpfsl_fb_wm  = NaN(size(fb_values));

%% ---------------------- Simulations (gray matter–like) ------------------
for i = 1:numel(T1a_values)
    T1a_curr = T1a_values(i);
    Rmpfsl_T1a_gm(i) = compute_Rmpfsl( ...
        T1a_curr, T1a_curr, T2a_gm, kba_gm, T2b_gm, fb_gm, ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

for i = 1:numel(T2a_values)
    Rmpfsl_T2a_gm(i) = compute_Rmpfsl( ...
        T1a_gm, T1a_gm, T2a_values(i), kba_gm, T2b_gm, fb_gm, ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

for i = 1:numel(kba_values)
    Rmpfsl_kba_gm(i) = compute_Rmpfsl( ...
        T1a_gm, T1a_gm, T2a_gm, kba_values(i), T2b_gm, fb_gm, ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

for i = 1:numel(T2b_values)
    Rmpfsl_T2b_gm(i) = compute_Rmpfsl( ...
        T1a_gm, T1a_gm, T2a_gm, kba_gm, T2b_values(i), fb_gm, ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

for i = 1:numel(fb_values)
    Rmpfsl_fb_gm(i) = compute_Rmpfsl( ...
        T1a_gm, T1a_gm, T2a_gm, kba_gm, T2b_gm, fb_values(i), ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

%% ---------------------- Simulations (white matter–like) -----------------
for i = 1:numel(T1a_values)
    T1a_curr = T1a_values(i);
    Rmpfsl_T1a_wm(i) = compute_Rmpfsl( ...
        T1a_curr, T1a_curr, T2a_wm, kba_wm, T2b_wm, fb_wm, ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

for i = 1:numel(T2a_values)
    Rmpfsl_T2a_wm(i) = compute_Rmpfsl( ...
        T1a_wm, T1a_wm, T2a_values(i), kba_wm, T2b_wm, fb_wm, ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

for i = 1:numel(kba_values)
    Rmpfsl_kba_wm(i) = compute_Rmpfsl( ...
        T1a_wm, T1a_wm, T2a_wm, kba_values(i), T2b_wm, fb_wm, ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

for i = 1:numel(T2b_values)
    Rmpfsl_T2b_wm(i) = compute_Rmpfsl( ...
        T1a_wm, T1a_wm, T2a_wm, kba_wm, T2b_values(i), fb_wm, ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

for i = 1:numel(fb_values)
    Rmpfsl_fb_wm(i) = compute_Rmpfsl( ...
        T1a_wm, T1a_wm, T2a_wm, kba_wm, T2b_wm, fb_values(i), ...
        Np, tp, td, w1_hi_scaled, dw_hi_scaled, w1_lo_scaled, dw_lo_scaled, ...
        Gz, tg, phi, angle, v, initSigns);
end

%% ---------------------- Plotting ----------------------------------------
figure('Color', 'w', ...
       'Name', 'R_{mpfsps} vs Various Parameters', ...
       'Units', 'pixels', ...
       'Position', [100 100 1100 480]);

% 使用 'none' 进一步减小间距
tl = tiledlayout(2, 3, 'TileSpacing', 'tight', 'Padding', 'tight');

% 使用指定的颜色：紫色和绿色
color_gm = [0.4940, 0.1840, 0.5560];   % 紫色 (Gray Matter)
color_wm = [0.4660, 0.6740, 0.1880];   % 绿色 (White Matter)

lineWidth = 1.8;
labelFontSize = 11;
titleFontSize = 11;
axisFontSize = 9;

% Plot T1a
ax = nexttile; hold(ax, 'on'); box(ax, 'on');
plot(ax, T1a_values * 1e3, Rmpfsl_T1a_gm, '-', 'LineWidth', lineWidth, 'Color', color_gm);
plot(ax, T1a_values * 1e3, Rmpfsl_T1a_wm, '-', 'LineWidth', lineWidth, 'Color', color_wm);
hold(ax, 'off');
xlabel(ax, 'T_{1a} (ms)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
ylabel(ax, 'R_{mpfsps} (s^{-1})', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title(ax, '(A) T_{1a}', 'FontSize', titleFontSize, 'FontWeight', 'bold');
ylim(ax, [0, 7]);
ax.FontSize = axisFontSize;
grid(ax, 'on');
legend(ax, {'Gray Matter', 'White Matter'}, ...
       'Location', 'northwest', 'Box', 'off', 'FontSize', 8);

% Plot T2a
ax = nexttile; hold(ax, 'on'); box(ax, 'on');
plot(ax, T2a_values * 1e3, Rmpfsl_T2a_gm, '-', 'LineWidth', lineWidth, 'Color', color_gm);
plot(ax, T2a_values * 1e3, Rmpfsl_T2a_wm, '-', 'LineWidth', lineWidth, 'Color', color_wm);
hold(ax, 'off');
xlabel(ax, 'T_{2a} (ms)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
% ylabel(ax, 'R_{mpfsps} (s^{-1})', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title(ax, '(B) T_{2a}', 'FontSize', titleFontSize, 'FontWeight', 'bold');
ylim(ax, [0, 8]);
ax.FontSize = axisFontSize;
grid(ax, 'on');

% Plot kba
ax = nexttile; hold(ax, 'on'); box(ax, 'on');
plot(ax, kba_values, Rmpfsl_kba_gm, '-', 'LineWidth', lineWidth, 'Color', color_gm);
plot(ax, kba_values, Rmpfsl_kba_wm, '-', 'LineWidth', lineWidth, 'Color', color_wm);
hold(ax, 'off');
xlabel(ax, 'k_{ba} (s^{-1})', 'FontSize', labelFontSize, 'FontWeight', 'bold');
% ylabel(ax, 'R_{mpfsps} (s^{-1})', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title(ax, '(C) k_{ba}', 'FontSize', titleFontSize, 'FontWeight', 'bold');
ylim(ax, [0, 13]);
ax.FontSize = axisFontSize;
grid(ax, 'on');

% Plot T2b
ax = nexttile; hold(ax, 'on'); box(ax, 'on');
plot(ax, T2b_values * 1e6, Rmpfsl_T2b_gm, '-', 'LineWidth', lineWidth, 'Color', color_gm);
plot(ax, T2b_values * 1e6, Rmpfsl_T2b_wm, '-', 'LineWidth', lineWidth, 'Color', color_wm);
hold(ax, 'off');
xlabel(ax, 'T_{2b} (\mus)', 'FontSize', labelFontSize, 'FontWeight', 'bold');
ylabel(ax, 'R_{mpfsps} (s^{-1})', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title(ax, '(D) T_{2b}', 'FontSize', titleFontSize, 'FontWeight', 'bold');
ylim(ax, [0, 7]);
ax.FontSize = axisFontSize;
grid(ax, 'on');

% Plot fb
ax = nexttile; hold(ax, 'on'); box(ax, 'on');
plot(ax, fb_values, Rmpfsl_fb_gm, '-', 'LineWidth', lineWidth, 'Color', color_gm);
plot(ax, fb_values, Rmpfsl_fb_wm, '-', 'LineWidth', lineWidth, 'Color', color_wm);
hold(ax, 'off');
xlabel(ax, 'f_b', 'FontSize', labelFontSize, 'FontWeight', 'bold');
% ylabel(ax, 'R_{mpfsps} (s^{-1})', 'FontSize', labelFontSize, 'FontWeight', 'bold');
title(ax, '(E) f_b', 'FontSize', titleFontSize, 'FontWeight', 'bold');
ylim(ax, [0, 10]);
ax.FontSize = axisFontSize;
grid(ax, 'on');

% Empty sixth tile - 隐藏
nexttile; axis off;

exportgraphics(gcf, 'figure water.png', 'Resolution', 600);

%% ========================================================================
%% ---------------------- Helper functions --------------------------------
%% ========================================================================
function Rmpfsl = compute_Rmpfsl(T1a, T1b, T2a, kba, T2b, fb, ...
                                 Np, tp, td, ...
                                 w1_hi, dw_hi, w1_lo, dw_lo, ...
                                 Gz, tg, phi, angle, v, initSigns)
    % Enforce R1b to track R1a by using the provided T1b (default = T1a).
    R1a = 1 / T1a;
    R2a = 1 / T2a;
    R1b = 1 / T1b;
    kab = fb * kba;
    M0a = 1;
    M0b = fb;

    [Mz_hi_pos, Mz_hi_neg, Mz_lo_pos, Mz_lo_neg] = ...
        simulate_four_experiments_epg(Np, tp, td, ...
                                      w1_hi, dw_hi, w1_lo, dw_lo, ...
                                      Gz, tg, phi, angle, ...
                                      R1a, R2a, R1b, ...
                                      kab, kba, M0a, M0b, T2b, v, initSigns);

    diff_hi = Mz_hi_pos - Mz_hi_neg;
    diff_lo = Mz_lo_pos - Mz_lo_neg;

    if diff_hi > 0 && diff_lo > 0
        Rmpfsl = log(diff_lo / diff_hi) / (Np * tp);
    else
        Rmpfsl = NaN;
    end
end

function [Mz_hi_pos, Mz_hi_neg, Mz_lo_pos, Mz_lo_neg] = ...
        simulate_four_experiments_epg(Np, tp, td, ...
                                      w1_hi, dw_hi, w1_lo, dw_lo, ...
                                      Gz, tg, phi, angle, ...
                                      R1a, R2a, R1b, ...
                                      kab, kba, M0a, M0b, T2b, v, initSigns)

    Mz_hi_pos = run_epg_final(Np, tp, td, w1_hi, +dw_hi, ...
                              Gz, tg, phi, angle, ...
                              R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v, initSigns(1));

    Mz_hi_neg = run_epg_final(Np, tp, td, w1_hi, -dw_hi, ...
                              Gz, tg, phi, angle, ...
                              R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v, initSigns(2));

    Mz_lo_pos = run_epg_final(Np, tp, td, w1_lo, +dw_lo, ...
                              Gz, tg, phi, angle, ...
                              R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v, initSigns(3));

    Mz_lo_neg = run_epg_final(Np, tp, td, w1_lo, -dw_lo, ...
                              Gz, tg, phi, angle, ...
                              R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v, initSigns(4));
end

function Mz_final = run_epg_final(Np, tp, td, w1, dw, ...
                                  Gz, tg, phi, angle, ...
                                  R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v, initSign)
    [theta, Omega_eff] = compute_rf_angles(w1, dw);
    Rrfb = RF_MT(T2b, abs(w1), abs(dw), 'SuperLorentzian');
    beta = Omega_eff * tp;

    S = [0; 0; initSign * M0a; initSign * M0b];

    for n = 1:Np
        S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
        S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
        if abs(Gz) > 0
            S = epgmt_Grad(S);
        end
        S = epgmt_Flow(S, Gz, tg, v, angle);
    end

    Mz_final = real(S(3));
end

function [theta, Omega_eff] = compute_rf_angles(w1, dw)
    Omega_eff = hypot(w1, dw);
    if Omega_eff < 1e-12
        theta = 0;
        Omega_eff = 0;
    else
        cos_theta = dw / Omega_eff;
        cos_theta = max(min(cos_theta, 1), -1);
        theta = acos(cos_theta);
    end
end