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
% fb_component controls which tissue uses the swept f_b values.
scenarios = {
    struct('primary', gm, 'primary_name', 'Gray Matter', 'primary_frac', 1.0, 'primary_v', 0.0, ...
           'blood_frac', 0.0, 'blood_v', 0.0, 'label', 'Pure Gray Matter', ...
           'fb_component', 'primary')
    struct('primary', gm, 'primary_name', 'Gray Matter', 'primary_frac', 0.7, 'primary_v', 0.0, ...
           'blood_frac', 0.3, 'blood_v', 0.0, 'label', '70% Gray Matter + 30% Blood', ...
           'fb_component', 'primary')
    struct('primary', gm, 'primary_name', 'Gray Matter', 'primary_frac', 0.7, 'primary_v', 0.0, ...
           'blood_frac', 0.3, 'blood_v', 10.0, 'label', '70% Gray Matter + 30% Blood (blood suppression)', ...
           'fb_component', 'primary')
    struct('primary', wm, 'primary_name', 'White Matter', 'primary_frac', 1.0, 'primary_v', 0.0, ...
           'blood_frac', 0.0, 'blood_v', 0.0, 'label', 'Pure White Matter', ...
           'fb_component', 'primary')
    struct('primary', wm, 'primary_name', 'White Matter', 'primary_frac', 0.7, 'primary_v', 0.0, ...
           'blood_frac', 0.3, 'blood_v', 0.0, 'label', '70% White Matter + 30% Blood', ...
           'fb_component', 'primary')
    struct('primary', wm, 'primary_name', 'White Matter', 'primary_frac', 0.7, 'primary_v', 0.0, ...
           'blood_frac', 0.3, 'blood_v', 10.0, 'label', '70% White Matter + 30% Blood (blood suppression)', ...
           'fb_component', 'primary')
};

n_scenarios = numel(scenarios);

%% Storage for EPG results
results_epg = NaN(n_scenarios, n_fb);

% Storage for combined Mz curves: dimensions (scenario, experiment, n_fb)
% Experiment order: [Hi+, Hi-, Lo+, Lo-]
Mz_curves = NaN(n_scenarios, 4, n_fb);

%% Labels for the four experiments
exp_labels = { ...
    'High power (+\omega_d)', ...
    'High power (-\omega_d)', ...
    'Low power (+\omega_d)',  ...
    'Low power (-\omega_d)'};

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
        
        % Store the combined Mz curves for later plotting
        Mz_curves(s, 1, f_idx) = Mz_hi_pos;
        Mz_curves(s, 2, f_idx) = Mz_hi_neg;
        Mz_curves(s, 3, f_idx) = Mz_lo_pos;
        Mz_curves(s, 4, f_idx) = Mz_lo_neg;
        
        % Compute R_mpfsl for the current f_b
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

%% Plotting R_mpfsl (EPG only) versus f_b - Combined Figure with PVE
is_gm = cellfun(@(sc) strcmpi(sc.primary_name, 'Gray Matter'), scenarios);
is_wm = cellfun(@(sc) strcmpi(sc.primary_name, 'White Matter'), scenarios);

colors = lines(n_scenarios);
line_styles = {'-', '--', ':', '-.', '-', '--', ':', '-.'};

% Create a single figure with 1x3 tiledlayout - 增加宽度以在右侧留白
figure('Color', 'w', 'Position', [100, 100, 1280, 400]);
t = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'tight');

%% Tile 1: Grey matter scenarios
nexttile;
hold on; grid on;
gm_idx = find(is_gm);
for k = 1:numel(gm_idx)
    s = gm_idx(k);
    ls = line_styles{mod(k-1, numel(line_styles))+1};
    plot(fb_values, results_epg(s, :), ...
        'LineWidth', 2, ...
        'Color', colors(k, :), ...
        'LineStyle', ls, ...
        'DisplayName', scenarios{s}.label);
end
title('(A) Gray Matter', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('f_b', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('R_{mpfsl} (s^{-1})', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'Box', 'off', 'FontSize', 9);
set(gca, 'FontSize', 10, 'LineWidth', 1);
box on;

%% Tile 2: White matter scenarios
nexttile;
hold on; grid on;
wm_idx = find(is_wm);
for k = 1:numel(wm_idx)
    s = wm_idx(k);
    ls = line_styles{mod(k-1, numel(line_styles))+1};
    plot(fb_values, results_epg(s, :), ...
        'LineWidth', 2, ...
        'Color', colors(k, :), ...
        'LineStyle', ls, ...
        'DisplayName', scenarios{s}.label);
end
title('(B) White Matter', 'FontWeight', 'bold', 'FontSize', 12);
xlabel('f_b', 'FontSize', 11, 'FontWeight', 'bold');
% ylabel('R_{mpfsl} (s^{-1})', 'FontSize', 11, 'FontWeight', 'bold');
legend('Location', 'best', 'Box', 'off', 'FontSize', 9);
set(gca, 'FontSize', 10, 'LineWidth', 1);
box on;

%% Tile 3: PVE Schematic
ax_pve = nexttile;
hold(ax_pve, 'on');

% 正方形边长为1，面积为1
% 红色区域需要占30%，即面积 = 0.3
% 使用直角三角形: 面积 = 0.5 * base * height = 0.3
% 如果 base = height = a, 则 0.5 * a^2 = 0.3, a = sqrt(0.6) ≈ 0.775
a = sqrt(0.6);  % ≈ 0.775

% 70% 灰质/白质区域 (蓝灰色)
gm_x = [0, 1, 1, 1-a, 0];
gm_y = [0, 0, 1-a, 1, 1];
patch(ax_pve, gm_x, gm_y, [0.70, 0.70, 0.82], 'EdgeColor', 'none');

% 30% 血液/CSF区域 (浅红色)
blood_x = [1-a, 1, 1, 1-a];
blood_y = [1, 1-a, 1, 1];
patch(ax_pve, blood_x, blood_y, [1, 0.82, 0.82], 'EdgeColor', 'none');

% 添加整齐的斜线纹理到血液区域
num_lines = 6;
for i = 1:num_lines
    t_line = i / (num_lines + 1);
    
    % 斜线平行于斜边
    x1 = (1-a) + t_line * a;
    y1 = 1;
    x2 = 1;
    y2 = (1-a) + t_line * a;
    
    plot(ax_pve, [x1, x2], [y1, y2], 'Color', [0.75, 0.2, 0.2], 'LineWidth', 1.2);
end

% 绘制边框
plot(ax_pve, [0,1,1,0,0], [0,0,1,1,0], 'k-', 'LineWidth', 2.5);

% 绘制分界线（斜边）
plot(ax_pve, [1-a, 1], [1, 1-a], 'k-', 'LineWidth', 1.5);

% 添加文字标签 - GM/WM
text(ax_pve, 0.3, 0.5, {'70%', 'GM/WM'}, 'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Color', [0.25, 0.25, 0.45]);

% 添加文字标签 - Blood/CSF 直接写在红色区域内
% 三角形重心位置
cx = (1-a + 1 + 1) / 3;
cy = (1 + 1-a + 1) / 3;

text(ax_pve, cx, cy, {'30%', 'Blood/', 'CSF'}, 'FontSize', 16, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'Color', [0.6, 0, 0], ...
    'VerticalAlignment', 'middle');

% 添加体素标签
text(ax_pve, 0.5, -0.12, 'Single Voxel', 'FontSize', 12, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');

% 设置坐标轴
axis(ax_pve, 'equal');
xlim(ax_pve, [-0.15, 1.15]);
ylim(ax_pve, [-0.2, 1.15]);
axis(ax_pve, 'off');

% 添加PVE标题
title(ax_pve, '(C) Partial Volume Effect', 'FontSize', 12, 'FontWeight', 'bold');

% 使用exportgraphics导出（推荐方法）
exportgraphics(gcf, 'figure pve.png', 'Resolution', 600);

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