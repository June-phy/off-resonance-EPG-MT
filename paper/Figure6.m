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
tp_default = 0.7e-3;           % RF pulse duration [s]
td_default = 2.0e-3;           % spacing between DANTE modules [s]
Np_curves  = 100;              % number of RF pulses
Gz         = 25e-3;            % gradient amplitude [T/m]
tg         = 1e-3;             % gradient duration [s]
phi        = 0;                % RF phase [rad]
angle      = 0;                % flow direction relative to gradient

%% ---------------------- Tissue parameters (3T) --------------------------
% Blood
tissue(1).name = 'Blood';
tissue(1).T1   = 1932e-3;      % [s]
tissue(1).T2   = 275e-3;       % [s]
tissue(1).T2b  = 280e-6;        % [s]
tissue(1).fb   = 0.028;         % macromolecular pool fraction
tissue(1).kba  = 35;           % exchange rate [s^-1]

% White Matter (WM)
tissue(2).name = 'WM';
tissue(2).T1   = 1084e-3;      % [s]
tissue(2).T2   = 69e-3;        % [s]
tissue(2).T2b  = 10e-6;        % [s]
tissue(2).fb   = 0.139;        % macromolecular pool fraction
tissue(2).kba  = 23;           % exchange rate [s^-1]

% Gray Matter (GM)
tissue(3).name = 'GM';
tissue(3).T1   = 1820e-3;      % [s]
tissue(3).T2   = 99e-3;        % [s]
tissue(3).T2b  = 9.1e-6;        % [s]
tissue(3).fb   = 0.095;        % macromolecular pool fraction
tissue(3).kba  = 40;           % exchange rate [s^-1]

n_tissue = numel(tissue);

% 计算派生参数
for i = 1:n_tissue
    tissue(i).R1a = 1 / tissue(i).T1;
    tissue(i).R2a = 1 / tissue(i).T2;
    tissue(i).R1b = tissue(i).R1a;
    tissue(i).kab = tissue(i).fb * tissue(i).kba;
    tissue(i).M0a = 1;
    tissue(i).M0b = tissue(i).fb;
end

%% ---------------------- Velocity parameter ------------------------------
v = 2e-2;  % 固定速度 2 cm/s

%% ---------------------- RF parameters (两组) ----------------------------
rf_params(1).w1 = 2*pi*100;    % w1 = 100 Hz
rf_params(1).dw = 2*pi*1000;   % dw = 1000 Hz
rf_params(1).label = 'w1=100Hz, dw=1000Hz';

rf_params(2).w1 = 2*pi*500;    % w1 = 500 Hz
rf_params(2).dw = 2*pi*5000;   % dw = 5000 Hz
rf_params(2).label = 'w1=500Hz, dw=5000Hz';

n_rf = numel(rf_params);

%% ---------------------- 颜色和线型设置 ----------------------------------
% 组织对应的颜色
tissueColors = [0.8500 0.3250 0.0980;   % 橙色 - Blood
                0 0.4470 0.7410;         % 蓝色 - WM
                0.4660 0.6740 0.1880];   % 绿色 - GM

% RF参数对应的线型
rfLineStyles = {'-', '--'};  % 第一组实线，第二组虚线
rfLineWidths = [2.5, 2.5];

%% ---------------------- Compute Mz curves -------------------------------
N_axis = 0:Np_curves;
Mz_v0 = zeros(n_tissue, n_rf, Np_curves + 1);
Mz_v  = zeros(n_tissue, n_rf, Np_curves + 1);
Ratio = zeros(n_tissue, n_rf, Np_curves + 1);

eps_ratio = 1e-9;

for rf_idx = 1:n_rf
    w1_rad = rf_params(rf_idx).w1;
    dw_rad = rf_params(rf_idx).dw;
    
    for t_idx = 1:n_tissue
        R1a = tissue(t_idx).R1a;
        R2a = tissue(t_idx).R2a;
        R1b = tissue(t_idx).R1b;
        kab = tissue(t_idx).kab;
        kba = tissue(t_idx).kba;
        M0a = tissue(t_idx).M0a;
        M0b = tissue(t_idx).M0b;
        T2b = tissue(t_idx).T2b;
        
        [theta, Omega_eff] = compute_rf_angles(w1_rad, dw_rad);
        Rrfb = RF_MT(T2b, w1_rad, dw_rad, 'SuperLorentzian');
        beta = Omega_eff * tp_default;
        
        % v=0 的曲线
        S = [0; 0; M0a; M0b];
        Mz_v0(t_idx, rf_idx, 1) = M0a;
        
        for n = 1:Np_curves
            S = epgmt_RF(S, theta, beta, phi, Rrfb, tp_default);
            S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp_default, td_default, M0a, M0b);
            S = epgmt_Grad(S);
            S = epgmt_Flow(S, Gz, tg, 0, angle);
            Mz_v0(t_idx, rf_idx, n + 1) = real(S(3));
        end
        
        % v=2cm/s 的曲线
        S = [0; 0; M0a; M0b];
        Mz_v(t_idx, rf_idx, 1) = M0a;
        
        for n = 1:Np_curves
            S = epgmt_RF(S, theta, beta, phi, Rrfb, tp_default);
            S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp_default, td_default, M0a, M0b);
            S = epgmt_Grad(S);
            S = epgmt_Flow(S, Gz, tg, v, angle);
            Mz_v(t_idx, rf_idx, n + 1) = real(S(3));
        end
        
        Ratio(t_idx, rf_idx, :) = squeeze(Mz_v(t_idx, rf_idx, :)) ./ ...
                                   (squeeze(Mz_v0(t_idx, rf_idx, :)) + eps_ratio);
    end
end

%% ---------------------- Plot --------------------------------------------
figure('Color', 'w', 'Position', [200, 455, 449, 395]);
hold on; grid on;

% 为每组RF参数添加shade（Blood和WM之间）
shadeAlpha = 0.15;

for rf_idx = 1:n_rf
    % Blood (t_idx=1) 和 WM (t_idx=2) 之间的shade
    blood_curve = squeeze(Ratio(1, rf_idx, :));
    wm_curve = squeeze(Ratio(2, rf_idx, :));
    
    % 使用fill绘制shade区域
    x_fill = [N_axis, fliplr(N_axis)];
    y_fill = [blood_curve', fliplr(wm_curve')];
    
    % shade颜色：混合Blood和WM的颜色，或用灰色
    shadeColor = [0.5, 0.5, 0.5];  % 灰色
    
    fill(x_fill, y_fill, shadeColor, ...
        'FaceAlpha', shadeAlpha, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

% 画曲线（不加入legend）
for t_idx = 1:n_tissue
    for rf_idx = 1:n_rf
        ratio_curve = squeeze(Ratio(t_idx, rf_idx, :));
        plot(N_axis, ratio_curve, ...
            'Color', tissueColors(t_idx, :), ...
            'LineStyle', rfLineStyles{rf_idx}, ...
            'LineWidth', rfLineWidths(rf_idx), ...
            'HandleVisibility', 'off');
    end
end

% 创建dummy lines用于legend
h_dummy = gobjects(n_tissue + n_rf, 1);
dummy_labels = cell(n_tissue + n_rf, 1);

% 组织颜色图例（用实线表示）
for t_idx = 1:n_tissue
    h_dummy(t_idx) = plot(NaN, NaN, '-', ...
        'Color', tissueColors(t_idx, :), ...
        'LineWidth', 2.5);
    dummy_labels{t_idx} = tissue(t_idx).name;
end

% RF参数线型图例（用黑色表示）
for rf_idx = 1:n_rf
    h_dummy(n_tissue + rf_idx) = plot(NaN, NaN, rfLineStyles{rf_idx}, ...
        'Color', 'k', ...
        'LineWidth', 2);
    dummy_labels{n_tissue + rf_idx} = rf_params(rf_idx).label;
end

xlabel('N_p', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Ratio (M_z^{v=2} / M_z^{v=0})', 'FontSize', 14, 'FontWeight', 'bold');

xlim([0, Np_curves]);
all_ratio = Ratio(:);
y_min = floor(min(all_ratio) * 10) / 10;
ylim([y_min, 1.05]);

legend(h_dummy, dummy_labels, 'Location', 'northeast', 'FontSize', 11, 'Box', 'off');
set(gca, 'FontSize', 12);

exportgraphics(gcf, 'figure tissue.png', 'Resolution', 600);

%% ---------------------- Helper functions --------------------------------
function [theta, Omega_eff] = compute_rf_angles(w1, dw)
    Omega_eff = sqrt(w1.^2 + dw.^2);
    if Omega_eff < 1e-12
        theta = 0;
        Omega_eff = 0;
    else
        cos_theta = dw / Omega_eff;
        cos_theta = max(min(cos_theta, 1), -1);
        theta = acos(cos_theta);
    end
end