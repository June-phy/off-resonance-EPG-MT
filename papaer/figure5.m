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

%% Constants & shared settings
gamma  = 2*pi*42.577e6;  % rad/(s*T)
tp     = 0.7e-3;         % [s] RF pulse duration
td     = 2.0e-3;         % [s] spacing between DANTE modules
Np_fixed = 200;          % fixed number of pulses
fprintf('Number of DANTE pulses fixed at N_p = %d.\n', Np_fixed);

Gz     = 25e-3;          % [T/m] gradient amplitude
tg     = 1e-3;           % [s] gradient duration
angle  = 0;              % [rad] flow direction

%% RF settings for R_mpfsl calculation
% Low power experiment
w1_lo_Hz = 100;
dw_lo_Hz = 1000;
w1_lo = 2*pi*w1_lo_Hz;
dw_lo = 2*pi*dw_lo_Hz;

% High power experiment
w1_hi_Hz = 500;
dw_hi_Hz = 5000;
w1_hi = 2*pi*w1_hi_Hz;
dw_hi = 2*pi*dw_hi_Hz;

% Initial magnetization signs for four experiments
initSigns = [+1, -1, +1, -1];

%% Grey matter parameters (gm-like)
gm.label = 'Gray Matter';
gm.T1    = 1820e-3;
gm.T2    = 99e-3;
gm.T2b   = 9.1e-6;
gm.fb    = 0.05;
gm.kba   = 40;

%% White matter parameters (wm-like)
wm.label = 'White Matter';
wm.T1    = 1084e-3;
wm.T2    = 69e-3;
wm.T2b   = 10e-6;
wm.fb    = 0.139;
wm.kba   = 23;

%% Blood parameters
blood.label = 'Blood';
blood.T1    = 1932e-3;
blood.T2    = 275e-3;
blood.T2b   = 280e-6;
blood.fb    = 0.028;
blood.kba   = 35;

%% Sweep settings for f_b
fb_values = linspace(0.01, 0.20, 16);  % adjust as needed
n_fb      = numel(fb_values);

%% Simulation scenarios (grey & white matter mixes with blood)
mix_fracs = [0.9 0.8 0.7 0.6 0.5];
scenarios = {};

% Pure GM
scenarios{end+1} = struct('primary', gm, 'primary_name', 'Gray Matter', 'primary_frac', 1.0, 'primary_v', 0.0, ...
       'blood_frac', 0.0, 'blood_v', 0.0, 'label', 'Pure Gray Matter', ...
       'fb_component', 'primary');

% Pure WM
scenarios{end+1} = struct('primary', wm, 'primary_name', 'White Matter', 'primary_frac', 1.0, 'primary_v', 0.0, ...
       'blood_frac', 0.0, 'blood_v', 0.0, 'label', 'Pure White Matter', ...
       'fb_component', 'primary');

% GM + Blood (v=0 & v=10)
for f = mix_fracs
    scenarios{end+1} = struct('primary', gm, 'primary_name', 'Gray Matter', 'primary_frac', f, 'primary_v', 0.0, ...
           'blood_frac', 1-f, 'blood_v', 0.0, ...
           'label', sprintf('%.0f%% Gray Matter + %.0f%% Blood', f*100, (1-f)*100), ...
           'fb_component', 'primary');
    scenarios{end+1} = struct('primary', gm, 'primary_name', 'Gray Matter', 'primary_frac', f, 'primary_v', 0.0, ...
           'blood_frac', 1-f, 'blood_v', 10.0, ...
           'label', sprintf('%.0f%% Gray Matter + %.0f%% Blood (blood suppression)', f*100, (1-f)*100), ...
           'fb_component', 'primary');
end

% WM + Blood (v=0 & v=10)
for f = mix_fracs
    scenarios{end+1} = struct('primary', wm, 'primary_name', 'White Matter', 'primary_frac', f, 'primary_v', 0.0, ...
           'blood_frac', 1-f, 'blood_v', 0.0, ...
           'label', sprintf('%.0f%% White Matter + %.0f%% Blood', f*100, (1-f)*100), ...
           'fb_component', 'primary');
    scenarios{end+1} = struct('primary', wm, 'primary_name', 'White Matter', 'primary_frac', f, 'primary_v', 0.0, ...
           'blood_frac', 1-f, 'blood_v', 10.0, ...
           'label', sprintf('%.0f%% White Matter + %.0f%% Blood (blood suppression)', f*100, (1-f)*100), ...
           'fb_component', 'primary');
end

n_scenarios = numel(scenarios);

%% Storage for EPG results
results_epg = NaN(n_scenarios, n_fb);
Mz_curves = NaN(n_scenarios, 4, n_fb);

%% EPG simulations for all scenarios across f_b
for s = 1:n_scenarios
    scenario = scenarios{s};
    fprintf('\n=== Simulating (EPG): %s ===\n', scenario.label);
    
    primary_template = scenario.primary;
    primary_frac     = scenario.primary_frac;
    primary_v        = scenario.primary_v;
    
    blood_frac = scenario.blood_frac;
    blood_v    = scenario.blood_v;
    
    tic;
    for f_idx = 1:n_fb
        fb_val = fb_values(f_idx);
        
        primary_params = primary_template;
        blood_params   = blood;
        
        switch lower(scenario.fb_component)
            case 'primary'
                primary_params.fb = fb_val;
            case 'blood'
                blood_params.fb = fb_val;
            case 'both'
                primary_params.fb = fb_val;
                blood_params.fb   = fb_val;
            otherwise
                error('Unknown fb_component option: %s', scenario.fb_component);
        end
        
        if primary_frac > 0
            [Mz_hi_pos_primary, Mz_hi_neg_primary, Mz_lo_pos_primary, Mz_lo_neg_primary] = ...
                simulate_four_exp_epg(primary_params, primary_v, Np_fixed, tp, td, ...
                w1_hi, dw_hi, w1_lo, dw_lo, gamma, Gz, tg, angle, initSigns);
        else
            Mz_hi_pos_primary = 0; Mz_hi_neg_primary = 0;
            Mz_lo_pos_primary = 0; Mz_lo_neg_primary = 0;
        end
        
        if blood_frac > 0
            [Mz_hi_pos_B, Mz_hi_neg_B, Mz_lo_pos_B, Mz_lo_neg_B] = ...
                simulate_four_exp_epg(blood_params, blood_v, Np_fixed, tp, td, ...
                w1_hi, dw_hi, w1_lo, dw_lo, gamma, Gz, tg, angle, initSigns);
        else
            Mz_hi_pos_B = 0; Mz_hi_neg_B = 0;
            Mz_lo_pos_B = 0; Mz_lo_neg_B = 0;
        end
        
        % Weighted combination
        Mz_hi_pos = primary_frac * Mz_hi_pos_primary + blood_frac * Mz_hi_pos_B;
        Mz_hi_neg = primary_frac * Mz_hi_neg_primary + blood_frac * Mz_hi_neg_B;
        Mz_lo_pos = primary_frac * Mz_lo_pos_primary + blood_frac * Mz_lo_pos_B;
        Mz_lo_neg = primary_frac * Mz_lo_neg_primary + blood_frac * Mz_lo_neg_B;
        
        Mz_curves(s, 1, f_idx) = Mz_hi_pos;
        Mz_curves(s, 2, f_idx) = Mz_hi_neg;
        Mz_curves(s, 3, f_idx) = Mz_lo_pos;
        Mz_curves(s, 4, f_idx) = Mz_lo_neg;
        
        diff_hi = Mz_hi_pos - Mz_hi_neg;
        diff_lo = Mz_lo_pos - Mz_lo_neg;
        
        if diff_hi > 0 && diff_lo > 0
            results_epg(s, f_idx) = log(diff_lo / diff_hi) / (Np_fixed * tp);
        else
            results_epg(s, f_idx) = NaN;
        end
    end
    
    fprintf('EPG simulation time: %.4f s\n', toc);
end

%% Precompute arrays for filtering
is_gm = cellfun(@(sc) strcmpi(sc.primary_name, 'Gray Matter'), scenarios);
is_wm = cellfun(@(sc) strcmpi(sc.primary_name, 'White Matter'), scenarios);

primary_frac_arr = cellfun(@(sc) sc.primary_frac, scenarios);
blood_frac_arr   = cellfun(@(sc) sc.blood_frac, scenarios);
blood_v_arr      = cellfun(@(sc) sc.blood_v, scenarios);

tol = 1e-6;
idx_pure_gm = find(is_gm & abs(primary_frac_arr-1)<tol & abs(blood_frac_arr)<tol, 1);
idx_pure_wm = find(is_wm & abs(primary_frac_arr-1)<tol & abs(blood_frac_arr)<tol, 1);

%% ---- Calibration curves (move up here for later bias plots) ----
R_pure_gm = results_epg(idx_pure_gm, :);
R_pure_wm = results_epg(idx_pure_wm, :);

[R_pure_gm_sorted, sort_idx_gm] = sort(R_pure_gm);
fb_for_gm_interp = fb_values(sort_idx_gm);
[R_pure_wm_sorted, sort_idx_wm] = sort(R_pure_wm);
fb_for_wm_interp = fb_values(sort_idx_wm);

[R_pure_gm_sorted, ia_gm] = unique(R_pure_gm_sorted, 'stable');
fb_for_gm_interp = fb_for_gm_interp(ia_gm);
[R_pure_wm_sorted, ia_wm] = unique(R_pure_wm_sorted, 'stable');
fb_for_wm_interp = fb_for_wm_interp(ia_wm);

%% Plotting settings
fs_title = 12; fs_label = 11; fs_axis = 10;
lw_main = 2; lw_zero = 1.5;
ratio_colors = lines(numel(mix_fracs));

%% ========================================================================
%% MERGED FIGURE: Figure1 + Figure3 in 2x3 layout with shared legends
%% ========================================================================
figure('Color', 'w', 'Position', [100, 100, 1300, 700]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

%% (A) PVE schematic
ax_pve = nexttile(1);
hold(ax_pve, 'on');

a = sqrt(0.6);
gm_x = [0, 1, 1, 1-a, 0];
gm_y = [0, 0, 1-a, 1, 1];
patch(ax_pve, gm_x, gm_y, [0.70, 0.70, 0.82], 'EdgeColor', 'none');

blood_x = [1-a, 1, 1, 1-a];
blood_y = [1, 1-a, 1, 1];
patch(ax_pve, blood_x, blood_y, [1, 0.82, 0.82], 'EdgeColor', 'none');

num_lines = 6;
for i = 1:num_lines
    t_line = i / (num_lines + 1);
    x1 = (1-a) + t_line * a; y1 = 1;
    x2 = 1; y2 = (1-a) + t_line * a;
    plot(ax_pve, [x1, x2], [y1, y2], 'Color', [0.75, 0.2, 0.2], 'LineWidth', 1.2);
end

plot(ax_pve, [0,1,1,0,0], [0,0,1,1,0], 'k-', 'LineWidth', 2.5);
plot(ax_pve, [1-a, 1], [1, 1-a], 'k-', 'LineWidth', 1.5);

text(ax_pve, 0.32, 0.45, {'Tissue', '(GM/WM)'}, 'FontSize', 13, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Color', [0.25, 0.25, 0.45]);

cx = (1-a + 1 + 1) / 3;
cy = (1 + 1-a + 1) / 3;
text(ax_pve, cx, cy+0.05, {'Blood/CSF'}, 'FontSize', 13, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Color', [0.6, 0, 0], 'VerticalAlignment', 'middle');

text(ax_pve, 0.5, -0.12, 'Single Voxel', 'FontSize', 11, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

axis(ax_pve, 'equal');
xlim(ax_pve, [-0.15, 1.15]);
ylim(ax_pve, [-0.2, 1.15]);
axis(ax_pve, 'off');
title(ax_pve, '(A) Partial Volume Effect', 'FontSize', fs_title, 'FontWeight', 'bold');

%% (B) GM R_mpfsl
ax_gm = nexttile(2);
hold(ax_gm, 'on'); grid(ax_gm, 'on');

if ~isempty(idx_pure_gm)
    plot(ax_gm, fb_values, results_epg(idx_pure_gm, :), 'LineWidth', lw_main, 'Color', [0 0 0]);
end

for i = 1:numel(mix_fracs)
    f = mix_fracs(i);
    c = ratio_colors(i, :);
    idx0  = find(is_gm & abs(primary_frac_arr-f)<tol & blood_v_arr==0  & blood_frac_arr>0, 1);
    idx10 = find(is_gm & abs(primary_frac_arr-f)<tol & blood_v_arr==10 & blood_frac_arr>0, 1);
    if ~isempty(idx0),  plot(ax_gm, fb_values, results_epg(idx0, :),  'LineWidth', lw_main, 'Color', c, 'LineStyle', '-'); end
    if ~isempty(idx10), plot(ax_gm, fb_values, results_epg(idx10, :), 'LineWidth', lw_main, 'Color', c, 'LineStyle', '--'); end
end

title(ax_gm, '(B) Gray Matter', 'FontWeight', 'bold', 'FontSize', fs_title);
xlabel(ax_gm, 'f_b', 'FontSize', fs_label, 'FontWeight', 'bold');
ylabel(ax_gm, 'R_{mpfsps} (s^{-1})', 'FontSize', fs_label, 'FontWeight', 'bold');
set(ax_gm, 'FontSize', fs_axis, 'LineWidth', 1); box(ax_gm, 'on');

%% (C) WM R_mpfsl
ax_wm = nexttile(3);
hold(ax_wm, 'on'); grid(ax_wm, 'on');

if ~isempty(idx_pure_wm)
    plot(ax_wm, fb_values, results_epg(idx_pure_wm, :), 'LineWidth', lw_main, 'Color', [0 0 0]);
end

for i = 1:numel(mix_fracs)
    f = mix_fracs(i);
    c = ratio_colors(i, :);
    idx0  = find(is_wm & abs(primary_frac_arr-f)<tol & blood_v_arr==0  & blood_frac_arr>0, 1);
    idx10 = find(is_wm & abs(primary_frac_arr-f)<tol & blood_v_arr==10 & blood_frac_arr>0, 1);
    if ~isempty(idx0),  plot(ax_wm, fb_values, results_epg(idx0, :),  'LineWidth', lw_main, 'Color', c, 'LineStyle', '-'); end
    if ~isempty(idx10), plot(ax_wm, fb_values, results_epg(idx10, :), 'LineWidth', lw_main, 'Color', c, 'LineStyle', '--'); end
end

title(ax_wm, '(C) White Matter', 'FontWeight', 'bold', 'FontSize', fs_title);
xlabel(ax_wm, 'f_b', 'FontSize', fs_label, 'FontWeight', 'bold');
set(ax_wm, 'FontSize', fs_axis, 'LineWidth', 1); box(ax_wm, 'on');

%% Legend tile (English)
ax_leg = nexttile(4);
hold(ax_leg, 'on'); axis(ax_leg, 'off');

leg_handles = [];
leg_labels  = {};

leg_handles(end+1) = plot(ax_leg, nan, nan, 'k-', 'LineWidth', lw_main);
leg_labels{end+1} = 'Pure tissue (GM/WM)';

for i = 1:numel(mix_fracs)
    f = mix_fracs(i);
    c = ratio_colors(i, :);
    leg_handles(end+1) = plot(ax_leg, nan, nan, '-', 'Color', c, 'LineWidth', lw_main);
    leg_labels{end+1} = sprintf('GM/WM %.0f%% + Blood %.0f%%', f*100, (1-f)*100);
end

leg_handles(end+1) = plot(ax_leg, nan, nan, '-',  'Color', [0.3 0.3 0.3], 'LineWidth', lw_main);
leg_labels{end+1}  = 'Solid: No flow suppression';

leg_handles(end+1) = plot(ax_leg, nan, nan, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', lw_main);
leg_labels{end+1}  = 'Dashed: Flow suppression';

leg_handles(end+1) = plot(ax_leg, nan, nan, 'k:', 'LineWidth', lw_zero);
leg_labels{end+1}  = 'Bias = 0';

lgd = legend(ax_leg, leg_handles, leg_labels, 'Location', 'northwest', ...
    'Box', 'on', 'FontSize', 11);
title(ax_leg, 'Legend', 'FontSize', fs_title, 'FontWeight', 'bold');

% --- Move legend closer to center of the legend tile ---
ax_leg.Units = 'normalized';
lgd.Units = 'normalized';
axpos = ax_leg.Position;
lgd.Position = [axpos(1) + 0.08*axpos(3), ...
                axpos(2) + 0.20*axpos(4), ...
                0.84*axpos(3), ...
                0.70*axpos(4)];

%% (D) GM bias with pair-wise shaded bands
ax_gm_bias = nexttile(5);
hold(ax_gm_bias, 'on'); grid(ax_gm_bias, 'on');
yline(ax_gm_bias, 0, 'k:', 'LineWidth', lw_zero);

for i = 1:numel(mix_fracs)
    f = mix_fracs(i);
    c = ratio_colors(i, :);
    idx0  = find(is_gm & abs(primary_frac_arr-f)<tol & blood_v_arr==0  & blood_frac_arr>0, 1);
    idx10 = find(is_gm & abs(primary_frac_arr-f)<tol & blood_v_arr==10 & blood_frac_arr>0, 1);
    if ~isempty(idx0) && ~isempty(idx10)
        fb_app0 = interp1(R_pure_gm_sorted, fb_for_gm_interp, results_epg(idx0, :), 'linear', 'extrap');
        fb_app10 = interp1(R_pure_gm_sorted, fb_for_gm_interp, results_epg(idx10, :), 'linear', 'extrap');
        fb_app0 = min(max(fb_app0, min(fb_values)), max(fb_values));
        fb_app10 = min(max(fb_app10, min(fb_values)), max(fb_values));
        bias0 = fb_app0 - fb_values;
        bias10 = fb_app10 - fb_values;
        valid = isfinite(bias0) & isfinite(bias10);
        if any(valid)
            shade = 1 - (1 - c) * 0.2; % lighter version of line color
            fh = fill(ax_gm_bias, [fb_values(valid), fliplr(fb_values(valid))], ...
                [bias0(valid), fliplr(bias10(valid))], shade, ...
                'FaceAlpha', 0.25, 'EdgeColor', 'none');
            uistack(fh, 'bottom');
        end
    end
end

for i = 1:numel(mix_fracs)
    f = mix_fracs(i);
    c = ratio_colors(i, :);
    idx0  = find(is_gm & abs(primary_frac_arr-f)<tol & blood_v_arr==0  & blood_frac_arr>0, 1);
    idx10 = find(is_gm & abs(primary_frac_arr-f)<tol & blood_v_arr==10 & blood_frac_arr>0, 1);
    if ~isempty(idx0)
        fb_app = interp1(R_pure_gm_sorted, fb_for_gm_interp, results_epg(idx0, :), 'linear', 'extrap');
        fb_app = min(max(fb_app, min(fb_values)), max(fb_values));
        plot(ax_gm_bias, fb_values, fb_app - fb_values, '-', 'Color', c, 'LineWidth', lw_main);
    end
    if ~isempty(idx10)
        fb_app = interp1(R_pure_gm_sorted, fb_for_gm_interp, results_epg(idx10, :), 'linear', 'extrap');
        fb_app = min(max(fb_app, min(fb_values)), max(fb_values));
        plot(ax_gm_bias, fb_values, fb_app - fb_values, '--', 'Color', c, 'LineWidth', lw_main);
    end
end

title(ax_gm_bias, '(D) Gray Matter', 'FontWeight', 'bold', 'FontSize', fs_title);
xlabel(ax_gm_bias, 'True f_b', 'FontSize', fs_label, 'FontWeight', 'bold');
ylabel(ax_gm_bias, 'Bias (Estimated f_b - True f_b)', 'FontSize', fs_label, 'FontWeight', 'bold');
set(ax_gm_bias, 'FontSize', fs_axis, 'LineWidth', 1); box(ax_gm_bias, 'on');

%% (E) WM bias with pair-wise shaded bands
ax_wm_bias = nexttile(6);
hold(ax_wm_bias, 'on'); grid(ax_wm_bias, 'on');
yline(ax_wm_bias, 0, 'k:', 'LineWidth', lw_zero);

for i = 1:numel(mix_fracs)
    f = mix_fracs(i);
    c = ratio_colors(i, :);
    idx0  = find(is_wm & abs(primary_frac_arr-f)<tol & blood_v_arr==0  & blood_frac_arr>0, 1);
    idx10 = find(is_wm & abs(primary_frac_arr-f)<tol & blood_v_arr==10 & blood_frac_arr>0, 1);
    if ~isempty(idx0) && ~isempty(idx10)
        fb_app0 = interp1(R_pure_wm_sorted, fb_for_wm_interp, results_epg(idx0, :), 'linear', 'extrap');
        fb_app10 = interp1(R_pure_wm_sorted, fb_for_wm_interp, results_epg(idx10, :), 'linear', 'extrap');
        fb_app0 = min(max(fb_app0, min(fb_values)), max(fb_values));
        fb_app10 = min(max(fb_app10, min(fb_values)), max(fb_values));
        bias0 = fb_app0 - fb_values;
        bias10 = fb_app10 - fb_values;
        valid = isfinite(bias0) & isfinite(bias10);
        if any(valid)
            shade = 1 - (1 - c) * 0.2;
            fh = fill(ax_wm_bias, [fb_values(valid), fliplr(fb_values(valid))], ...
                [bias0(valid), fliplr(bias10(valid))], shade, ...
                'FaceAlpha', 0.25, 'EdgeColor', 'none');
            uistack(fh, 'bottom');
        end
    end
end

for i = 1:numel(mix_fracs)
    f = mix_fracs(i);
    c = ratio_colors(i, :);
    idx0  = find(is_wm & abs(primary_frac_arr-f)<tol & blood_v_arr==0  & blood_frac_arr>0, 1);
    idx10 = find(is_wm & abs(primary_frac_arr-f)<tol & blood_v_arr==10 & blood_frac_arr>0, 1);
    if ~isempty(idx0)
        fb_app = interp1(R_pure_wm_sorted, fb_for_wm_interp, results_epg(idx0, :), 'linear', 'extrap');
        fb_app = min(max(fb_app, min(fb_values)), max(fb_values));
        plot(ax_wm_bias, fb_values, fb_app - fb_values, '-', 'Color', c, 'LineWidth', lw_main);
    end
    if ~isempty(idx10)
        fb_app = interp1(R_pure_wm_sorted, fb_for_wm_interp, results_epg(idx10, :), 'linear', 'extrap');
        fb_app = min(max(fb_app, min(fb_values)), max(fb_values));
        plot(ax_wm_bias, fb_values, fb_app - fb_values, '--', 'Color', c, 'LineWidth', lw_main);
    end
end

title(ax_wm_bias, '(E) White Matter', 'FontWeight', 'bold', 'FontSize', fs_title);
xlabel(ax_wm_bias, 'True f_b', 'FontSize', fs_label, 'FontWeight', 'bold');
% ylabel(ax_wm_bias, 'Bias (Estimated - True)', 'FontSize', fs_label, 'FontWeight', 'bold');
set(ax_wm_bias, 'FontSize', fs_axis, 'LineWidth', 1); box(ax_wm_bias, 'on');

exportgraphics(gcf, 'figure pve revision.png', 'Resolution', 600);

%% ========================================================================
%% Helper Functions - EPG
%% ========================================================================

function [Mz_hi_pos, Mz_hi_neg, Mz_lo_pos, Mz_lo_neg] = ...
    simulate_four_exp_epg(params, v, Np, tp, td, ...
    w1_hi, dw_hi, w1_lo, dw_lo, gamma, Gz, tg, angle, initSigns)

    Mz_hi_pos = run_epg_exp(params, v, Np, tp, td, w1_hi, +dw_hi, ...
        gamma, Gz, tg, angle, initSigns(1));
    
    Mz_hi_neg = run_epg_exp(params, v, Np, tp, td, w1_hi, -dw_hi, ...
        gamma, Gz, tg, angle, initSigns(2));
    
    Mz_lo_pos = run_epg_exp(params, v, Np, tp, td, w1_lo, +dw_lo, ...
        gamma, Gz, tg, angle, initSigns(3));
    
    Mz_lo_neg = run_epg_exp(params, v, Np, tp, td, w1_lo, -dw_lo, ...
        gamma, Gz, tg, angle, initSigns(4));
end

function Mz_final = run_epg_exp(params, v, Np, tp, td, w1, dw, gamma, Gz, tg, angle, initSign)
    T1  = params.T1;
    T2  = params.T2;
    T2b = params.T2b;
    fb  = params.fb;
    kba = params.kba;
    kab = fb * kba;
    
    R1a = 1/T1;
    R2a = 1/T2;
    R1b = R1a;
    
    M0a = 1;
    M0b = fb;
    
    Rrfb = RF_MT(T2b, abs(w1), abs(dw), 'SuperLorentzian');
    [theta, Omega_eff] = compute_rf_angles(w1, dw);
    beta = Omega_eff * tp;
    
    S = [0; 0; initSign * M0a; initSign * M0b];
    
    for n = 1:Np
        S = epgmt_RF(S, theta, beta, 0, Rrfb, tp);
        S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
        S = epgmt_Grad(S);
        S = epgmt_Flow(S, Gz, tg, v, angle);
    end
    
    Mz_final = real(S(3, 1));
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