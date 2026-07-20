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

%% ---------------------- Common parameters -------------------------------
gamma      = 2.6752e8;
tp_default = 0.7e-3;
td_default = 2.0e-3;
Np_curves  = 350;
Np_maps    = 100;
Np_used    = 100;
Gz         = 25e-3;
tg         = 1e-3;
phi        = 0;
angle      = 0;

% Blood @ 3T
T1   = 1932e-3;
T2   = 275e-3;
T2b  = 280e-6;
fb   = 0.028;
kba  = 35;
kab  = fb * kba;
R1a  = 1 / T1;
R2a  = 1 / T2;
R1b  = R1a;
M0a  = 1;
M0b  = fb;

% Velocities
v_arr_cm = [0, 0.1, 0.5, 1.0, 5.0, 10.0, 40.0, 80.0];
v_arr    = v_arr_cm * 1e-2;
n_v      = numel(v_arr);

% -------- Custom velocity color palette (fixed: purple / blue / teal / green / orange / red)
v_colors =  [0.49 0.18 0.56;   % purple
            0.00 0.45 0.74;   % blue
            [255, 182, 193] / 255;
            0.47 0.67 0.19;   % green
            0.93 0.69 0.13;   % orange
            0.85 0.10 0.10;   % red
            [255, 247, 0] / 255; 
             [0.6902, 0.7686, 0.8706]];  


% -------- Operating points -----------------------------------------------
op_dw_Hz   = [1000, 5000];
op_w1_Hz   = [ 100,  500];
op_labels_dwW1 = { ...
    {'\omega_1 = 100 Hz','\Delta\omega = 1000 Hz'}, ...
    {'\omega_1 = 500 Hz','\Delta\omega = 5000 Hz'} };

op_tp_ms = 0.7;
op_td_ms = 2.0;
op_labels_tpTd = { {'t_p = 0.7 ms','t_d = 2 ms'} };

%% ---------------------- TSE readout settings ----------------------------
TSE.N_startup       = 3;
TSE.ETL_acq         = 72;
TSE.ETL_tot         = TSE.N_startup + TSE.ETL_acq;
TSE.ESP             = 3.4e-3;
TSE.tp_ex           = 0.95e-3;
TSE.tp_ref          = 0.78e-3;
TSE.flip_ex         = pi/2;
TSE.flip_ref        = pi;
TSE.phi_ex          = 0;
TSE.phi_ref         = pi/2;
TSE.theta_tse       = pi/2;
TSE.dw_tse          = 0;
TSE.w1_ex           = TSE.flip_ex  / TSE.tp_ex;
TSE.w1_ref          = TSE.flip_ref / TSE.tp_ref;
TSE.half_esp        = TSE.ESP / 2;
TSE.gap_ex_to_ref1  = TSE.ESP / 2;
TSE.echo_idx_TEeff  = TSE.N_startup + 1;
TSE.TEeff_ms        = TSE.echo_idx_TEeff * TSE.ESP * 1e3;

TSE.flow_enable     = true;
TSE.Gz_readout      = Gz;
TSE.tg_readout      = tg;

TSE.Rrfb_ex  = RF_MT(T2b, TSE.w1_ex,  TSE.dw_tse, 'SuperLorentzian');
TSE.Rrfb_ref = RF_MT(T2b, TSE.w1_ref, TSE.dw_tse, 'SuperLorentzian');

PAR.R1a=R1a; PAR.R2a=R2a; PAR.R1b=R1b;
PAR.kab=kab; PAR.kba=kba;
PAR.M0a=M0a; PAR.M0b=M0b;

fprintf('TSE readout: ETL=%d (+%d startup), ESP=%.2f ms, TE_eff=%.2f ms, flow=%s\n', ...
        TSE.ETL_acq, TSE.N_startup, TSE.ESP*1e3, TSE.TEeff_ms, string(TSE.flow_enable));

%% ---------------------- 1) Signal curves (w1=100 Hz) --------------------
fprintf('Subplot A: signal vs Np (w1=100 Hz) ...\n'); tic;
w1_case1 = 2*pi*100;
dw_case1 = 2*pi*1000;
Sig_curves_case1 = simulate_velocity_curves(w1_case1, dw_case1, tp_default, td_default, ...
                                            Np_curves, Gz, tg, phi, angle, ...
                                            PAR, T2b, v_arr, TSE);
toc;

%% ---------------------- 2) Signal map vs (Δω, ω1) at v=0 ----------------
fprintf('Subplot B: Mxy map (dw,w1) v=0 ...\n'); tic;
dw_Hz_vals = linspace(0, 15000, 121);
w1_Hz_vals = linspace(0, 1000, 101);
Sig_dw_w1_v0 = compute_dw_w1_map(dw_Hz_vals, w1_Hz_vals, tp_default, td_default, ...
                                 Np_maps, Gz, tg, phi, angle, ...
                                 PAR, T2b, 0, TSE);
toc;

%% ---------------------- 3) Signal map vs (t_p, t_d) at v=0 --------------
fprintf('Subplot E: Mxy map (tp,td) v=0 ...\n'); tic;
tp_ms = linspace(0.01, 2.0, 80);
td_ms = linspace(0.01, 10.0, 100);
w1_base = 2*pi*100;
dw_base = 2*pi*1000;
[theta_base, Omega_eff_base] = compute_rf_angles(w1_base, dw_base);
Rrfb_base = RF_MT(T2b, w1_base, dw_base, 'SuperLorentzian');
Sig_tp_td_v0 = compute_tp_td_map(tp_ms, td_ms, Np_maps, Gz, tg, phi, angle, ...
                                 PAR, theta_base, Omega_eff_base, Rrfb_base, 0, TSE);
toc;

%% ---------------------- 4) Signal curves (w1=500 Hz) --------------------
fprintf('Subplot D: signal vs Np (w1=500 Hz) ...\n'); tic;
w1_case4 = 2*pi*500;
dw_case4 = 2*pi*5000;
Sig_curves_case4 = simulate_velocity_curves(w1_case4, dw_case4, tp_default, td_default, ...
                                            Np_curves, Gz, tg, phi, angle, ...
                                            PAR, T2b, v_arr, TSE);
toc;

%% ---------------------- 5) v=1 signal map (dw,w1) & ratio --------------
fprintf('Subplot C: Mxy map (dw,w1) v=1 cm/s ...\n'); tic;
Sig_dw_w1_v10 = compute_dw_w1_map(dw_Hz_vals, w1_Hz_vals, tp_default, td_default, ...
                                  Np_maps, Gz, tg, phi, angle, ...
                                  PAR, T2b, 1e-2, TSE);
toc;
eps_ratio = 1e-9;
Sig_ratio_dw_w1 = Sig_dw_w1_v10 ./ (Sig_dw_w1_v0 + eps_ratio);

%% ---------------------- 6) v=1 signal map (tp,td) & ratio --------------
fprintf('Subplot F: Mxy map (tp,td) v=1 cm/s ...\n'); tic;
Sig_tp_td_v10 = compute_tp_td_map(tp_ms, td_ms, Np_maps, Gz, tg, phi, angle, ...
                                  PAR, theta_base, Omega_eff_base, Rrfb_base, 1e-2, TSE);
toc;
Sig_ratio_tp_td = Sig_tp_td_v10 ./ (Sig_tp_td_v0 + eps_ratio);

%% ---------------------- Assemble figure (2 × 3) -------------------------
N_axis_curves = 0:Np_curves;

bigFig = figure('Color', 'w', ...
                'Name', 'DANTE + TSE readout: Blood parameter selection (2x3)', ...
                'Position', [ 551         261        1013         640]);
tl = tiledlayout(bigFig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ylab_signal = 'M_{xy}';

color_lvl_0p1 = [0.10 0.85 0.10];
color_lvl_0p2 = [1.00 0.55 0.00];
color_lvl_0p3 = [1.00 0.20 0.20];

% Unified color limits [0,1] for all heatmaps
clim_unified = [0, 1];

contour_levels = [0.1, 0.2, 0.3];
contour_colors = [color_lvl_0p1; color_lvl_0p2; color_lvl_0p3];
contour_lw     = 0.9;

% --- (A) signal vs N (w1=100Hz) -----------------------------------------
ax1 = nexttile(tl, 1);
hold(ax1, 'on'); grid(ax1, 'on');
h_lines_A = gobjects(n_v, 1);
for j = 1:n_v
    h_lines_A(j) = plot(ax1, N_axis_curves, Sig_curves_case1(j, :), ...
                        'LineWidth', 1.6, 'Color', v_colors(j, :));
end
xl_A = xline(ax1, Np_used, '--', 'Color', [0.35 0.35 0.35], ...
             'LineWidth', 1.2, 'Label', sprintf('N_p = %d', Np_used), ...
             'LabelHorizontalAlignment', 'center', ...
             'LabelVerticalAlignment', 'top', ...
             'LabelOrientation', 'horizontal', ...
             'FontSize', 8, 'Interpreter', 'tex');
xl_A.HandleVisibility = 'off';
leg1 = legend(ax1, h_lines_A, ...
              arrayfun(@(vc) sprintf('v = %.1f cm/s', vc), v_arr_cm, ...
                       'UniformOutput', false), ...
              'Location', 'northeast', 'Box', 'off');
leg1.FontSize = 8;
xlabel(ax1, 'N_p', 'FontSize', 10, 'Interpreter', 'tex');
ylabel(ax1, ylab_signal, 'FontSize', 10);
title(ax1, '(A) \omega_1=100 Hz, \Delta\omega=1000 Hz', ...
      'FontSize', 10, 'FontWeight', 'bold');
xlim(ax1, [0, Np_curves]);
ax1.FontSize = 9;

% --- (B) Mxy heatmap (dw,w1) v=0 --------------------------------------
ax2 = nexttile(tl, 2);
imagesc(ax2, dw_Hz_vals, w1_Hz_vals, Sig_dw_w1_v0);
axis(ax2, 'tight');
set(ax2, 'YDir', 'normal', 'FontSize', 9);
colormap(ax2, parula); clim(ax2, clim_unified);
mark_operating_points(ax2, op_dw_Hz, op_w1_Hz, op_labels_dwW1);
xlabel(ax2, '\Delta\omega (Hz)', 'FontSize', 10);
text(ax2, -0.02, 1.02, '\omega_1 (Hz)', 'Units','normalized', ...
     'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom', ...
     'Interpreter','tex');
title(ax2, '(B) M_{xy}  (v = 0 cm/s)', 'FontSize', 10, 'FontWeight','bold');

% --- (C) Ratio heatmap (dw,w1) v=1/v=0 + ratio contours -----------------
ax3 = nexttile(tl, 3);
imagesc(ax3, dw_Hz_vals, w1_Hz_vals, Sig_ratio_dw_w1);
axis(ax3, 'tight');
set(ax3, 'YDir', 'normal', 'FontSize', 9);
colormap(ax3, parula); clim(ax3, clim_unified);
hold(ax3, 'on');
h_contours_C = gobjects(numel(contour_levels), 1);
for k = 1:numel(contour_levels)
    [~, h_tmp] = contour(ax3, dw_Hz_vals, w1_Hz_vals, Sig_ratio_dw_w1, ...
                         [contour_levels(k), contour_levels(k)], ...
                         'LineColor', contour_colors(k, :), ...
                         'LineWidth', contour_lw);
    h_contours_C(k) = h_tmp;
end
mark_operating_points(ax3, op_dw_Hz, op_w1_Hz, op_labels_dwW1);
xlabel(ax3, '\Delta\omega (Hz)', 'FontSize', 10);
text(ax3, -0.02, 1.02, '\omega_1 (Hz)', 'Units','normalized', ...
     'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom', ...
     'Interpreter','tex');
title(ax3, '(C) Ratio  (v = 1 cm/s / v = 0 cm/s)', ...
      'FontSize', 10, 'FontWeight','bold');
legend(ax3, h_contours_C, arrayfun(@(lvl) sprintf('ratio = %.1f', lvl), ...
       contour_levels, 'UniformOutput', false), ...
       'Location', 'northwest', 'FontSize', 8, ...
       'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');

% --- (D) signal vs N (w1=500Hz) -----------------------------------------
ax4 = nexttile(tl, 4);
hold(ax4, 'on'); grid(ax4, 'on');
for j = 1:n_v
    plot(ax4, N_axis_curves, Sig_curves_case4(j, :), ...
         'LineWidth', 1.6, 'Color', v_colors(j, :));
end
xl_D = xline(ax4, Np_used, '--', 'Color', [0.35 0.35 0.35], ...
             'LineWidth', 1.2, 'Label', sprintf('N_p = %d', Np_used), ...
             'LabelHorizontalAlignment', 'center', ...
             'LabelVerticalAlignment', 'top', ...
             'LabelOrientation', 'horizontal', ...
             'FontSize', 8, 'Interpreter', 'tex');
xl_D.HandleVisibility = 'off';
xlabel(ax4, 'N_p', 'FontSize', 10, 'Interpreter', 'tex');
ylabel(ax4, ylab_signal, 'FontSize', 10);
title(ax4, '(D) \omega_1=500 Hz, \Delta\omega=5000 Hz', ...
      'FontSize', 10, 'FontWeight', 'bold');
xlim(ax4, [0, Np_curves]);
ax4.FontSize = 9;

% --- (E) Mxy heatmap (tp,td) v=0 ----------------------------------------
ax5 = nexttile(tl, 5);
imagesc(ax5, tp_ms, td_ms, Sig_tp_td_v0);
axis(ax5, 'tight');
set(ax5, 'YDir', 'normal', 'FontSize', 9);
colormap(ax5, parula); clim(ax5, clim_unified);
mark_operating_points(ax5, op_tp_ms, op_td_ms, op_labels_tpTd);
xlabel(ax5, 't_p (ms)', 'FontSize', 10);
text(ax5, -0.02, 1.02, 't_d (ms)', 'Units','normalized', ...
     'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom', ...
     'Interpreter','tex');
title(ax5, '(E) M_{xy}  (v = 0 cm/s)', 'FontSize', 10, 'FontWeight','bold');

% --- (F) Ratio heatmap (tp,td) v=1/v=0 + ratio contours -----------------
ax6 = nexttile(tl, 6);
imagesc(ax6, tp_ms, td_ms, Sig_ratio_tp_td);
axis(ax6, 'tight');
set(ax6, 'YDir', 'normal', 'FontSize', 9);
colormap(ax6, parula); clim(ax6, clim_unified);
hold(ax6, 'on');
h_contours_F = gobjects(numel(contour_levels), 1);
for k = 1:numel(contour_levels)
    [~, h_tmp] = contour(ax6, tp_ms, td_ms, Sig_ratio_tp_td, ...
                         [contour_levels(k), contour_levels(k)], ...
                         'LineColor', contour_colors(k, :), ...
                         'LineWidth', contour_lw);
    h_contours_F(k) = h_tmp;
end
mark_operating_points(ax6, op_tp_ms, op_td_ms, op_labels_tpTd);
xlabel(ax6, 't_p (ms)', 'FontSize', 10);
text(ax6, -0.02, 1.02, 't_d (ms)', 'Units','normalized', ...
     'FontSize',10,'HorizontalAlignment','right','VerticalAlignment','bottom', ...
     'Interpreter','tex');
title(ax6, '(F) Ratio  (v = 1 cm/s / v = 0 cm/s)', ...
      'FontSize', 10, 'FontWeight','bold');
legend(ax6, h_contours_F, arrayfun(@(lvl) sprintf('ratio = %.1f', lvl), ...
       contour_levels, 'UniformOutput', false), ...
       'Location', 'northeast', 'FontSize', 8, ...
       'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');

% --- Shared colorbar attached to tiledlayout (east side) -----------------
cb = colorbar(ax6, 'eastoutside');
cb.Layout.Tile = 'east';
cb.FontSize    = 10;
cb.Ticks       = 0:0.2:1;
cb.Label.String   = 'M_{xy} / Ratio';
cb.Label.FontSize = 11;
cb.Label.Interpreter = 'tex';

exportgraphics(gcf, 'figure parameter selection2_revision.png', 'Resolution', 600);
fprintf('All subplots done.\n');


%% ========================================================================
%   Helper: mark operating points on heatmap using parameter-value labels
%% ========================================================================
function mark_operating_points(ax, xs, ys, labels)
    hold(ax, 'on');
    xl = xlim(ax); yl = ylim(ax);
    dx =  0.02 * (xl(2)-xl(1));
    dy =  0.03 * (yl(2)-yl(1));
    for k = 1:numel(xs)
        plot(ax, xs(k), ys(k), 'o', ...
             'MarkerSize', 9, ...
             'MarkerFaceColor', 'w', ...
             'MarkerEdgeColor', 'k', ...
             'LineWidth', 1.3, ...
             'HandleVisibility', 'off');
        plot(ax, xs(k), ys(k), '.', ...
             'MarkerSize', 9, ...
             'Color', 'k', ...
             'HandleVisibility', 'off');

        lab = labels{k};
        text(ax, xs(k)+dx, ys(k)+dy, lab, ...
             'Color', 'k', 'FontWeight', 'bold', 'FontSize', 9, ...
             'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
             'Interpreter','tex');
        text(ax, xs(k)+dx*0.85, ys(k)+dy*1.10, lab, ...
             'Color', 'w', 'FontWeight', 'bold', 'FontSize', 9, ...
             'HorizontalAlignment','left', 'VerticalAlignment','bottom', ...
             'Interpreter','tex');
    end
end

%% ========================================================================
%                            TSE readout core
%% ========================================================================
function sig = run_TSE_readout(S, PAR, TSE, v, flow_angle)
    if nargin < 4, v = 0;          end
    if nargin < 5, flow_angle = 0; end
    do_flow = TSE.flow_enable && (v ~= 0);

    S = epgmt_RF(S, TSE.theta_tse, TSE.flip_ex, TSE.phi_ex, ...
                 TSE.Rrfb_ex, TSE.tp_ex);
    S = epgmt_Ex_Relax(S, PAR.R1a, PAR.R2a, PAR.R1b, PAR.kab, PAR.kba, ...
                       0, TSE.gap_ex_to_ref1, PAR.M0a, PAR.M0b);
    if do_flow
        S = epgmt_Flow(S, TSE.Gz_readout, TSE.tg_readout, v, flow_angle);
    end

    sig = NaN;
    for e = 1:TSE.ETL_tot
        if e > 1
            S = epgmt_Ex_Relax(S, PAR.R1a, PAR.R2a, PAR.R1b, ...
                               PAR.kab, PAR.kba, 0, TSE.half_esp, ...
                               PAR.M0a, PAR.M0b);
            if do_flow
                S = epgmt_Flow(S, TSE.Gz_readout, TSE.tg_readout, v, flow_angle);
            end
        end
        S = epgmt_Grad(S);
        S = epgmt_RF(S, TSE.theta_tse, TSE.flip_ref, TSE.phi_ref, ...
                     TSE.Rrfb_ref, TSE.tp_ref);
        S = epgmt_Grad(S);
        S = epgmt_Ex_Relax(S, PAR.R1a, PAR.R2a, PAR.R1b, ...
                           PAR.kab, PAR.kba, 0, TSE.half_esp, ...
                           PAR.M0a, PAR.M0b);
        if do_flow
            S = epgmt_Flow(S, TSE.Gz_readout, TSE.tg_readout, v, flow_angle);
        end
        if e == TSE.echo_idx_TEeff
            sig = real(S(1,1));
            return;
        end
    end
end

%% ========================================================================
%                          Simulation helpers
%% ========================================================================
function Sig_curves = simulate_velocity_curves(w1, dw, tp, td, Np, Gz, tg, ...
                                               phi, angle, PAR, T2b, v_arr, TSE)
    n_v = numel(v_arr);
    Sig_curves = zeros(n_v, Np + 1);
    [theta, Omega_eff] = compute_rf_angles(w1, dw);
    Rrfb = RF_MT(T2b, w1, dw, 'SuperLorentzian');
    beta = Omega_eff * tp;

    for j = 1:n_v
        v = v_arr(j);
        S = [0; 0; PAR.M0a; PAR.M0b];
        Sig_curves(j, 1) = run_TSE_readout(S, PAR, TSE, v, angle);
        for n = 1:Np
            S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
            S = epgmt_Ex_Relax(S, PAR.R1a, PAR.R2a, PAR.R1b, ...
                               PAR.kab, PAR.kba, tp, td, PAR.M0a, PAR.M0b);
            S = epgmt_Grad(S);
            S = epgmt_Flow(S, Gz, tg, v, angle);
            Sig_curves(j, n + 1) = run_TSE_readout(S, PAR, TSE, v, angle);
        end
    end
end

% ------------------------------------------------------------------------
function Sig_map = compute_dw_w1_map(dw_Hz_vals, w1_Hz_vals, tp, td, Np, ...
                                     Gz, tg, phi, angle, PAR, T2b, v, TSE)
    n_dw = numel(dw_Hz_vals);
    n_w1 = numel(w1_Hz_vals);
    Sig_map = zeros(n_w1, n_dw);
    for iy = 1:n_w1
        w1_rad = 2*pi*w1_Hz_vals(iy);
        for ix = 1:n_dw
            dw_rad = 2*pi*dw_Hz_vals(ix);
            [theta, Omega_eff] = compute_rf_angles(w1_rad, dw_rad);
            Rrfb = RF_MT(T2b, w1_rad, dw_rad, 'SuperLorentzian');
            beta = Omega_eff * tp;
            S = [0; 0; PAR.M0a; PAR.M0b];
            for n = 1:Np
                S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
                S = epgmt_Ex_Relax(S, PAR.R1a, PAR.R2a, PAR.R1b, ...
                                   PAR.kab, PAR.kba, tp, td, PAR.M0a, PAR.M0b);
                S = epgmt_Grad(S);
                S = epgmt_Flow(S, Gz, tg, v, angle);
            end
            Sig_map(iy, ix) = run_TSE_readout(S, PAR, TSE, v, angle);
        end
        if mod(iy, 20) == 0
            fprintf('  dw_w1_map row %d / %d\n', iy, n_w1);
        end
    end
end

% ------------------------------------------------------------------------
function Sig_map = compute_tp_td_map(tp_ms, td_ms, Np, Gz, tg, phi, angle, ...
                                     PAR, theta_base, Omega_eff_base, ...
                                     Rrfb_base, v, TSE)
    n_tp = numel(tp_ms);
    n_td = numel(td_ms);
    Sig_map = zeros(n_td, n_tp);
    for iy = 1:n_td
        td = td_ms(iy) * 1e-3;
        for ix = 1:n_tp
            tp = tp_ms(ix) * 1e-3;
            beta = Omega_eff_base * tp;
            S = [0; 0; PAR.M0a; PAR.M0b];
            for n = 1:Np
                S = epgmt_RF(S, theta_base, beta, phi, Rrfb_base, tp);
                S = epgmt_Ex_Relax(S, PAR.R1a, PAR.R2a, PAR.R1b, ...
                                   PAR.kab, PAR.kba, tp, td, PAR.M0a, PAR.M0b);
                S = epgmt_Grad(S);
                S = epgmt_Flow(S, Gz, tg, v, angle);
            end
            Sig_map(iy, ix) = run_TSE_readout(S, PAR, TSE, v, angle);
        end
        if mod(iy, 20) == 0
            fprintf('  tp_td_map row %d / %d\n', iy, n_td);
        end
    end
end

% ------------------------------------------------------------------------
function [theta, Omega_eff] = compute_rf_angles(w1, dw)
    Omega_eff = sqrt(w1.^2 + dw.^2);
    if Omega_eff < 1e-12
        theta = 0;
        Omega_eff = 0;
    else
        cos_theta = max(min(dw / Omega_eff, 1), -1);
        theta = acos(cos_theta);
    end
end