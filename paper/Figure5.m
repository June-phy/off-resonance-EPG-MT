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

%% ---------------------- Common parameters -------------------------------
gamma      = 2.6752e8;         % rad/(s*T), kept for completeness
tp_default = 0.7e-3;           % RF pulse duration [s]
td_default = 2.0e-3;           % spacing between DANTE modules [s]
Np_curves  = 350;              % number of RF pulses for EPG curves (subplot A, D)
Np_maps   = 100;              % number of RF pulses for maps (subplot B, C, E, F)
Gz         = 25e-3;            % gradient amplitude [T/m]
tg         = 1e-3;             % gradient duration [s]
phi        = 0;                % RF phase [rad]
angle      = 0;                % flow direction relative to gradient

% Blood relaxation / MT parameters (3T)
T1   = 1932e-3;                % [s] - Blood T1 at 3T
T2   = 275e-3;                 % [s] - Blood T2 at 3T
T2b  = 280e-6;                  % [s] - macromolecular pool T2
fb   = 0.028;                   % macromolecular pool fraction (blood has less MT)
kba  = 35;                     % exchange rate [s^-1]
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
                                           Np_curves, Gz, tg, phi, angle, ...
                                           R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v_arr);

%% ---------------------- 2) Mz map vs (Δω, ω1) at v=0 --------------------
dw_Hz_vals = linspace(0, 15000, 121);   % Δω sweep [Hz]
w1_Hz_vals = linspace(0, 1000, 101);    % ω1 sweep [Hz]
Mz_dw_w1_v0 = compute_dw_w1_map(dw_Hz_vals, w1_Hz_vals, tp_default, td_default, ...
                                Np_maps, Gz, tg, phi, angle, ...
                                R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, 0);

%% ---------------------- 3) Mz map vs (t_p, t_d) at v=0 ------------------
tp_ms = linspace(0.01, 2.0, 80);       % RF duration sweep [ms]
td_ms = linspace(0.01, 10.0, 100);     % spacing sweep [ms]
w1_base = 2*pi*100;                     % rad/s (constant for this scan)
dw_base = 2*pi*1000;                    % rad/s
[theta_base, Omega_eff_base] = compute_rf_angles(w1_base, dw_base);
Rrfb_base = RF_MT(T2b, w1_base, dw_base, 'SuperLorentzian');
Mz_tp_td_v0 = compute_tp_td_map(tp_ms, td_ms, Np_maps, Gz, tg, phi, angle, ...
                                R1a, R2a, R1b, kab, kba, M0a, M0b, ...
                                theta_base, Omega_eff_base, Rrfb_base, 0);

%% ---------------------- 4) EPG curves (w1=500 Hz) -----------------------
w1_case4 = 2*pi*500;           % rad/s
dw_case4 = 2*pi*5000;          % rad/s
Mz_curves_case4 = simulate_velocity_curves(w1_case4, dw_case4, tp_default, td_default, ...
                                           Np_curves, Gz, tg, phi, angle, ...
                                           R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, v_arr);

%% ---------------------- 5) Ratio map vs (Δω, ω1) ------------------------
Mz_dw_w1_v10 = compute_dw_w1_map(dw_Hz_vals, w1_Hz_vals, tp_default, td_default, ...
                                 Np_maps, Gz, tg, phi, angle, ...
                                 R1a, R2a, R1b, kab, kba, M0a, M0b, T2b, 1e-2);
eps_ratio = 1e-9;
Mz_ratio_dw_w1 = Mz_dw_w1_v10 ./ (Mz_dw_w1_v0 + eps_ratio);

%% ---------------------- 6) Ratio map vs (t_p, t_d) ----------------------
Mz_tp_td_v10 = compute_tp_td_map(tp_ms, td_ms, Np_maps, Gz, tg, phi, angle, ...
                                 R1a, R2a, R1b, kab, kba, M0a, M0b, ...
                                 theta_base, Omega_eff_base, Rrfb_base, 1e-2);
Mz_ratio_tp_td = Mz_tp_td_v10 ./ (Mz_tp_td_v0 + eps_ratio);

%% ---------------------- Assemble figure (2 × 3 tiles) -------------------
N_axis_curves = 0:Np_curves;

bigFig = figure('Color', 'w', ...
                'Name', 'EPG DANTE MT Summary - Blood (2×3)', ...
                'Position', [655, 262, 1011, 627]);

% 使用 'tight' 进一步减小间距
tl = tiledlayout(bigFig, 2, 3, 'TileSpacing', 'tight', 'Padding', 'compact');

% --- Subplot 1 (A): Mz vs N (w1=100 Hz, Δω=1000 Hz) ---
ax1 = nexttile(tl, 1);
hold(ax1, 'on'); grid(ax1, 'on');
for j = 1:n_v
    plot(ax1, N_axis_curves, Mz_curves_case1(j, :), 'LineWidth', 1.5);
end
leg1 = legend(ax1, arrayfun(@(vc) sprintf('v = %.1f cm/s', vc), v_arr_cm, 'UniformOutput', false), ...
              'Location', 'southwest', 'Box', 'off');
leg1.FontSize = 8;
xlabel(ax1, 'N_p', 'FontSize', 10, 'Interpreter', 'tex');
ylabel(ax1, 'M_z', 'FontSize', 10);
title(ax1, '(A) EPG (\omega_1 = 100 Hz, \Delta\omega = 1000 Hz)', ...
      'FontSize', 10, 'FontWeight', 'bold');
xlim(ax1, [0, Np_curves]);
ax1.FontSize = 9;

% --- Subplot 2 (B): Mz heatmap (Δω, ω1) @ v=0 ---
ax2 = nexttile(tl, 2);
imagesc(ax2, dw_Hz_vals, w1_Hz_vals, Mz_dw_w1_v0);
axis(ax2, 'tight');
set(ax2, 'YDir', 'normal', 'FontSize', 9);
colormap(ax2, parula);
clim(ax2, [0, 1]);
xlabel(ax2, '\Delta\omega (Hz)', 'FontSize', 10);
% 将 ylabel 移到 y 轴顶部
ylabel(ax2, '');  % 清空默认 ylabel
text(ax2, -0.02, 1.02, '\omega_1 (Hz)', 'Units', 'normalized', ...
     'FontSize', 10, 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom', 'Interpreter', 'tex');
title(ax2, '(B) M_z (v = 0 cm/s)', ...
      'FontSize', 10, 'FontWeight', 'bold');

% --- Subplot 3 (C): Ratio heatmap (Δω, ω1) ---
ax3 = nexttile(tl, 3);
imagesc(ax3, dw_Hz_vals, w1_Hz_vals, Mz_ratio_dw_w1);
axis(ax3, 'tight');
set(ax3, 'YDir', 'normal', 'FontSize', 9);
colormap(ax3, parula);
clim(ax3, [0, 1]);
hold(ax3, 'on');
contour_levels_C = [0.3, 0.2, 0.1];
contour_colors_C = [0 1 1; 1 0 0; 0 0 1];
h_contours_C = gobjects(numel(contour_levels_C), 1);
for k = 1:numel(contour_levels_C)
    [~, h_tmp] = contour(ax3, dw_Hz_vals, w1_Hz_vals, Mz_ratio_dw_w1, ...
                         [contour_levels_C(k), contour_levels_C(k)], ...
                         'LineColor', contour_colors_C(k, :), ...
                         'LineWidth', 1.8);
    h_contours_C(k) = h_tmp;
end
xlabel(ax3, '\Delta\omega (Hz)', 'FontSize', 10);
% 将 ylabel 移到 y 轴顶部
ylabel(ax3, '');
text(ax3, -0.02, 1.02, '\omega_1 (Hz)', 'Units', 'normalized', ...
     'FontSize', 10, 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom', 'Interpreter', 'tex');
title(ax3, '(C) Ratio (v = 1 vs 0 cm/s)', ...
      'FontSize', 10, 'FontWeight', 'bold');
legend(ax3, h_contours_C, arrayfun(@(lvl) sprintf('%.1f', lvl), ...
       contour_levels_C, 'UniformOutput', false), ...
       'Location', 'northwest', 'FontSize', 8);

% --- Subplot 4 (D): Mz vs N (w1=500 Hz, Δω=5000 Hz) ---
ax4 = nexttile(tl, 4);
hold(ax4, 'on'); grid(ax4, 'on');
for j = 1:n_v
    plot(ax4, N_axis_curves, Mz_curves_case4(j, :), 'LineWidth', 1.5);
end
xlabel(ax4, 'N_p', 'FontSize', 10, 'Interpreter', 'tex');
ylabel(ax4, 'M_z', 'FontSize', 10);
title(ax4, '(D) EPG (\omega_1 = 500 Hz, \Delta\omega = 5000 Hz)', ...
      'FontSize', 10, 'FontWeight', 'bold');
xlim(ax4, [0, Np_curves]);
ax4.FontSize = 9;

% --- Subplot 5 (E): Mz heatmap (t_p, t_d) @ v=0 ---
ax5 = nexttile(tl, 5);
imagesc(ax5, tp_ms, td_ms, Mz_tp_td_v0);
axis(ax5, 'tight');
set(ax5, 'YDir', 'normal', 'FontSize', 9);
colormap(ax5, parula);
clim(ax5, [0, 1]);
xlabel(ax5, 't_p (ms)', 'FontSize', 10);
% 将 ylabel 移到 y 轴顶部
ylabel(ax5, '');
text(ax5, -0.02, 1.02, 't_d (ms)', 'Units', 'normalized', ...
     'FontSize', 10, 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom', 'Interpreter', 'tex');
title(ax5, '(E) M_z (v = 0 cm/s)', ...
      'FontSize', 10, 'FontWeight', 'bold');

% --- Subplot 6 (F): Ratio heatmap (t_p, t_d) ---
ax6 = nexttile(tl, 6);
imagesc(ax6, tp_ms, td_ms, Mz_ratio_tp_td);
axis(ax6, 'tight');
set(ax6, 'YDir', 'normal', 'FontSize', 9);
colormap(ax6, parula);
clim(ax6, [0, 1]);
hold(ax6, 'on');
contour_levels_F = [0.1, 0.2, 0.3];
contour_colors_F = [0.1 0.7 0.1; 1.0 0.6 0.0; 0.8 0.0 0.2];
h_contours_F = gobjects(numel(contour_levels_F), 1);
for k = 1:numel(contour_levels_F)
    [~, h_tmp] = contour(ax6, tp_ms, td_ms, Mz_ratio_tp_td, ...
                         [contour_levels_F(k), contour_levels_F(k)], ...
                         'LineColor', contour_colors_F(k, :), ...
                         'LineWidth', 1.8);
    h_contours_F(k) = h_tmp;
end
xlabel(ax6, 't_p (ms)', 'FontSize', 10);
% 将 ylabel 移到 y 轴顶部
ylabel(ax6, '');
text(ax6, -0.02, 1.02, 't_d (ms)', 'Units', 'normalized', ...
     'FontSize', 10, 'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'bottom', 'Interpreter', 'tex');
title(ax6, '(F) Ratio (v = 1 cm/s vs 0 cm/s)', ...
      'FontSize', 10, 'FontWeight', 'bold');
legend(ax6, h_contours_F, arrayfun(@(lvl) sprintf('%.1f', lvl), contour_levels_F, ...
       'UniformOutput', false), 'Location', 'northwest', 'FontSize', 8);

% --- 添加共享的colorbar在最右侧 ---
cb = colorbar(ax6, 'eastoutside');
cb.Layout.Tile = 'east';  % 放在tiledlayout的最右侧
cb.FontSize = 10;
cb.Label.String = 'M_z / Ratio';
cb.Label.FontSize = 11;
cb.Ticks = 0:0.2:1;

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

exportgraphics(gcf, 'figure parameter selection2.png', 'Resolution', 600);