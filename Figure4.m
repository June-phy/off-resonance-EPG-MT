%-------------------------------------------------------------------------
% Script:      EPG DANTE MT Simulation Summary (2×3 Figure Layout)
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
%   This script simulates the effects of magnetization transfer (MT) 
%   in liver-like tissues under DANTE pulse trains using the Extended 
%   Phase Graph (EPG) framework. Key features include:
%
%     - Simulation of EPG curves for different velocities.
%     - Computation of magnetization (Mz) maps as a function of off-resonance 
%       frequency (Δω) and RF nutation rate (ω₁) at v = 0 and v = 10 cm/s.
%     - Computation of Mz maps as a function of RF pulse duration (tp) and 
%       spacing (td) for v = 0 and v = 10 cm/s.
%     - Generation of a 2×3 tiled figure showing magnetization dynamics 
%       under different conditions.
%
%   Subplots:
%     (a) Mz curves vs Np (w1 = 100 Hz, Δω = 1000 Hz).
%     (b) Mz heatmap (Δω, ω₁) for v = 0 cm/s.
%     (c) Mz heatmap (tp, td) for v = 0 cm/s.
%     (d) Mz curves vs Np (w1 = 400 Hz, Δω = 4000 Hz).
%     (e) Ratio map (Δω, ω1) for v = 10 cm/s vs v = 0 cm/s.
%     (f) Ratio map (tp, td) for v = 10 cm/s vs v = 0 cm/s.
%
% Inputs:
%   - Relaxation and exchange parameters for liver tissue.
%   - Velocity array for flow simulation.
%   - RF pulse parameters (duration, spacing, nutation rates).
%
% Outputs:
%   - A 2×3 tiled figure displaying the simulated magnetization curves
%     and maps under various experimental conditions.
%-------------------------------------------------------------------------
%%
clear; close all; clc;

%% ---------------------- Common parameters -------------------------------
gamma      = 2.6752e8;         %#ok<NASGU> % rad/(s*T), kept for completeness
tp_default = 0.7e-3;           % RF pulse duration [s]
td_default = 1.0e-3;           % spacing between DANTE modules [s]
Np         = 200;              % number of RF pulses
Gz         = 5e-3;             % gradient amplitude [T/m]
tg         = 1e-3;             % gradient duration [s]
phi        = 0;                % RF phase [rad]
angle      = 0;                % flow direction relative to gradient

% Tissue (liver-like) relaxation / MT parameters
T1   = 812e-3;                 % [s]
T2   = 42e-3;                  % [s]
T2b  = 7.7e-6;                 % [s]
fb   = 0.069;                  % macromolecular pool fraction
kba  = 51;                     % exchange rate [s^-1]
kab  = fb * kba;

R1a  = 1 / T1;
R2a  = 1 / T2;
R1b  = R1a;

M0a  = 1;
M0b  = fb;

% Velocity array (cm/s -> m/s)
v_arr_cm = [0, 0.1, 0.5, 1.0, 5.0, 10.0];
v_arr    = v_arr_cm * 1e-2;    % m/s
n_v      = numel(v_arr);

%% ---------------------- 1) EPG curves (w1=100 Hz) -----------------------
w1_case1 = 2*pi*100;           % rad/s
dw_case1 = 2*pi*1000;          % rad/s
Mz_curves_case1 = simulate_velocity_curves(w1_case1, dw_case1, tp_default, td_default, ...
                                           Np, Gz, tg, phi, angle, ...
                                           R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v_arr);

%% ---------------------- 2) Mz map vs (Δω, ω1) at v=0 --------------------
dw_Hz_vals = linspace(0, 5000, 121);   % Δω sweep [Hz]
w1_Hz_vals = linspace(0, 1000, 101);   % ω1 sweep [Hz]
Mz_dw_w1_v0 = compute_dw_w1_map(dw_Hz_vals, w1_Hz_vals, tp_default, td_default, ...
                                Np, Gz, tg, phi, angle, ...
                                R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, 0);

%% ---------------------- 3) Mz map vs (t_p, t_d) at v=0 ------------------
tp_ms = linspace(0.01, 2.0, 80);       % RF duration sweep [ms]
td_ms = linspace(0.01, 10.0, 100);     % spacing sweep [ms]
w1_base = 2*pi*100;                     % rad/s (constant for this scan)
dw_base = 2*pi*1000;                    % rad/s
[theta_base, Omega_eff_base] = compute_rf_angles(w1_base, dw_base);
Rrfb_base = RF_MT(T2b, w1_base, dw_base, 'SuperLorentzian');
Mz_tp_td_v0 = compute_tp_td_map(tp_ms, td_ms, Np, Gz, tg, phi, angle, ...
                                R1a, R2a, R1b, kab, kba, M0a, M0b, ...
                                theta_base, Omega_eff_base, Rrfb_base, 0);

%% ---------------------- 4) EPG curves (w1=400 Hz) -----------------------
w1_case4 = 2*pi*400;           % rad/s
dw_case4 = 2*pi*4000;          % rad/s
Mz_curves_case4 = simulate_velocity_curves(w1_case4, dw_case4, tp_default, td_default, ...
                                           Np, Gz, tg, phi, angle, ...
                                           R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v_arr);

%% ---------------------- 5) Ratio map vs (Δω, ω1) ------------------------
Mz_dw_w1_v10 = compute_dw_w1_map(dw_Hz_vals, w1_Hz_vals, tp_default, td_default, ...
                                 Np, Gz, tg, phi, angle, ...
                                 R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, 10e-2);
eps_ratio = 1e-9;
Mz_ratio_dw_w1 = Mz_dw_w1_v10 ./ (Mz_dw_w1_v0 + eps_ratio);

%% ---------------------- 6) Ratio map vs (t_p, t_d) ----------------------
Mz_tp_td_v10 = compute_tp_td_map(tp_ms, td_ms, Np, Gz, tg, phi, angle, ...
                                 R1a, R2a, R1b, kab, kba, M0a, M0b, ...
                                 theta_base, Omega_eff_base, Rrfb_base, 10e-2);
Mz_ratio_tp_td = Mz_tp_td_v10 ./ (Mz_tp_td_v0 + eps_ratio);

%% ---------------------- Assemble figure (2 × 3 tiles) -------------------
N_axis = 0:Np;

bigFig = figure('Color', 'w', ...
                'Name', 'EPG DANTE MT Summary (2×3)', ...
                'Position', [100, 100, 1600, 900]);

tl = tiledlayout(bigFig, 2, 3, 'TileSpacing', 'tight', 'Padding', 'compact');

% --- Subplot 1: Mz vs N (w1=100 Hz, Δω=1000 Hz) ---
ax1 = nexttile(tl, 1);
hold(ax1, 'on'); grid(ax1, 'on');
for j = 1:n_v
    plot(ax1, N_axis, Mz_curves_case1(j, :), 'LineWidth', 1.5);
end
leg1 = legend(ax1, arrayfun(@(vc) sprintf('v = %.1f cm/s', vc), v_arr_cm, 'UniformOutput', false), ...
              'Location', 'southwest', 'Box', 'off');
leg1.FontSize = 11;
xlabel(ax1, 'N_p', 'FontSize', 12, 'Interpreter', 'tex');
ylabel(ax1, 'M_z', 'FontSize', 12);
title(ax1, '(a) EPG (w_1 = 100 Hz, \Delta\omega = 1000 Hz)', ...
      'FontSize', 13, 'FontWeight', 'bold');
xlim(ax1, [0, Np]);
ax1.FontSize = 12;

% --- Subplot 2: Mz heatmap (Δω, ω1) @ v=0 ---
ax2 = nexttile(tl, 2);
imagesc(ax2, dw_Hz_vals, w1_Hz_vals, Mz_dw_w1_v0);
axis(ax2, 'tight');
set(ax2, 'YDir', 'normal', 'FontSize', 12);
colormap(ax2, parula);
cb2 = colorbar(ax2);
cb2.Label.String = '';
xlabel(ax2, '\Delta\omega (Hz)', 'FontSize', 12);
ylabel(ax2, '\omega_1 (Hz)', 'FontSize', 12);
title(ax2, sprintf('(b) M_z (v = %.1f cm/s)', 0), ...
      'FontSize', 13, 'FontWeight', 'bold');

% --- Subplot 3: Mz heatmap (t_p, t_d) @ v=0 ---
ax3 = nexttile(tl, 3);
imagesc(ax3, tp_ms, td_ms, Mz_tp_td_v0);
axis(ax3, 'tight');
set(ax3, 'YDir', 'normal', 'FontSize', 12);
colormap(ax3, parula);
cb3 = colorbar(ax3);
cb3.Label.String = 'M_z';
cb3.Label.FontSize = 11;
xlabel(ax3, 't_p (ms)', 'FontSize', 12);
ylabel(ax3, 't_d (ms)', 'FontSize', 12);
title(ax3, sprintf('(c) M_z (v = %.1f cm/s)', 0), ...
      'FontSize', 13, 'FontWeight', 'bold');

% --- Subplot 4: Mz vs N (w1=400 Hz, Δω=4000 Hz) ---
ax4 = nexttile(tl, 4);
hold(ax4, 'on'); grid(ax4, 'on');
for j = 1:n_v
    plot(ax4, N_axis, Mz_curves_case4(j, :), 'LineWidth', 1.5);
end
xlabel(ax4, 'N_p', 'FontSize', 12, 'Interpreter', 'tex');
ylabel(ax4, 'M_z', 'FontSize', 12);
title(ax4, '(d) EPG (w_1 = 400 Hz, \Delta\omega = 4000 Hz)', ...
      'FontSize', 13, 'FontWeight', 'bold');
xlim(ax4, [0, Np]);
ax4.FontSize = 12;

% --- Subplot 5: Ratio heatmap (Δω, ω1) ---
ax5 = nexttile(tl, 5);
imagesc(ax5, dw_Hz_vals, w1_Hz_vals, Mz_ratio_dw_w1);
axis(ax5, 'tight');
set(ax5, 'YDir', 'normal', 'FontSize', 12);
colormap(ax5, parula);
cb5 = colorbar(ax5);
% cb5.Label.String = 'v = 10 cm/s / v = 0 cm/s';
cb5.Label.FontSize = 11;
hold(ax5, 'on');
[~, h_contour5] = contour(ax5, dw_Hz_vals, w1_Hz_vals, Mz_ratio_dw_w1, [0.3 0.3], ...
                          'LineColor', [0 1 1], 'LineWidth', 2.0);
xlabel(ax5, '\Delta\omega (Hz)', 'FontSize', 12);
ylabel(ax5, '\omega_1 (Hz)', 'FontSize', 12);
title(ax5, '(e) Ratio map (v = 10 cm/s vs 0 cm/s)', ...
      'FontSize', 13, 'FontWeight', 'bold');
leg5 = legend(ax5, h_contour5, 'ratio = 0.3', 'Location', 'northwest', 'FontSize', 11);
% leg5.TextColor = [0 1 1];
leg5.EdgeColor = 'none';

% --- Subplot 6: Ratio heatmap (t_p, t_d) ---
ax6 = nexttile(tl, 6);
imagesc(ax6, tp_ms, td_ms, Mz_ratio_tp_td);
axis(ax6, 'tight');
set(ax6, 'YDir', 'normal', 'FontSize', 12);
colormap(ax6, parula);
cb6 = colorbar(ax6);
cb6.Label.String = 'ratio';
cb6.Label.FontSize = 11;
hold(ax6, 'on');

contour_levels = [0.1, 0.2, 0.3];
contour_colors = [0.1 0.7 0.1;
                  1.0 0.6 0.0;
                  0.8 0.0 0.2];
h_contours = gobjects(numel(contour_levels), 1);

for k = 1:numel(contour_levels)
    [~, h_tmp] = contour(ax6, tp_ms, td_ms, Mz_ratio_tp_td, ...
                         [contour_levels(k), contour_levels(k)], ...
                         'LineColor', contour_colors(k, :), ...
                         'LineWidth', 1.8);
    h_contours(k) = h_tmp;
end

xlabel(ax6, 't_p (ms)', 'FontSize', 12);
ylabel(ax6, 't_d (ms)', 'FontSize', 12);
title(ax6, '(f) Ratio map (v = 10 cm/s vs 0 cm/s)', ...
      'FontSize', 13, 'FontWeight', 'bold');
legend(ax6, h_contours, arrayfun(@(lvl) sprintf('ratio = %.1f', lvl), contour_levels, ...
                                'UniformOutput', false), 'Location', 'northwest', 'FontSize', 11);

%% ---------------------- Helper functions --------------------------------
function Mz_curves = simulate_velocity_curves(w1, dw, tp, td, ...
                                              Np, Gz, tg, phi, angle, ...
                                              R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v_arr)
    n_v = numel(v_arr);
    Mz_curves = zeros(n_v, Np + 1);
    [theta, Omega_eff] = compute_rf_angles(w1, dw);
    Rrfb = RF_MT(T2b, w1, dw, 'SuperLorentzian');
    beta = Omega_eff * tp;

    for j = 1:n_v
        v = v_arr(j);
        S = [0; 0; M0a; M0b];
        Mz_curves(j, 1) = M0a;

        for n = 1:Np
            S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
            S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
            S = epgmt_Grad(S);
            S = epgmt_Flow(S, Gz, tg, v, angle);

            Mz_curves(j, n + 1) = real(S(3));
        end
    end
end

function Mz_map = compute_dw_w1_map(dw_Hz_vals, w1_Hz_vals, tp, td, ...
                                    Np, Gz, tg, phi, angle, ...
                                    R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v)
    n_dw = numel(dw_Hz_vals);
    n_w1 = numel(w1_Hz_vals);
    Mz_map = zeros(n_w1, n_dw);

    for iy = 1:n_w1
        w1_Hz = w1_Hz_vals(iy);
        w1_rad = 2*pi*w1_Hz;

        for ix = 1:n_dw
            dw_Hz = dw_Hz_vals(ix);
            dw_rad = 2*pi*dw_Hz;

            [theta, Omega_eff] = compute_rf_angles(w1_rad, dw_rad);
            Rrfb = RF_MT(T2b, w1_rad, dw_rad, 'SuperLorentzian');
            beta = Omega_eff * tp;

            S = [0; 0; M0a; M0b];

            for n = 1:Np
                S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
                S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
                S = epgmt_Grad(S);
                S = epgmt_Flow(S, Gz, tg, v, angle);
            end

            Mz_map(iy, ix) = real(S(3));
        end
    end
end

function Mz_map = compute_tp_td_map(tp_ms, td_ms, Np, Gz, tg, phi, angle, ...
                                    R1a, R2a, R1b, kab, kba, M0a, M0b, ...
                                    theta_base, Omega_eff_base, Rrfb_base, v)
    n_tp = numel(tp_ms);
    n_td = numel(td_ms);
    Mz_map = zeros(n_td, n_tp);

    for iy = 1:n_td
        td = td_ms(iy) * 1e-3;  % seconds

        for ix = 1:n_tp
            tp = tp_ms(ix) * 1e-3;  % seconds
            beta = Omega_eff_base * tp;

            S = [0; 0; M0a; M0b];

            for n = 1:Np
                S = epgmt_RF(S, theta_base, beta, phi, Rrfb_base, tp);
                S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
                S = epgmt_Grad(S);
                S = epgmt_Flow(S, Gz, tg, v, angle);
            end

            Mz_map(iy, ix) = real(S(3));
        end
    end
end

function [theta, Omega_eff] = compute_rf_angles(w1, dw)
    Omega_eff = sqrt(w1.^2 + dw.^2);
    if Omega_eff < 1e-12
        theta = 0;
        Omega_eff = 0;
    else
        cos_theta = dw / Omega_eff;
        cos_theta = max(min(cos_theta, 1), -1); % numerical safety
        theta = acos(cos_theta);
    end
end