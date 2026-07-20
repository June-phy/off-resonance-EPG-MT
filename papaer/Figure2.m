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

clear; close all;

%% 
gamma = 2.6752*1e8;

% 开关：控制字体是否加粗
use_bold_font = true;  % true = 加粗, false = 正常

if use_bold_font
    font_weight = 'bold';
else
    font_weight = 'normal';
end

T1a = 1084*1e-3; 
T2a = 69*1e-3;
T2b = 10*1.e-6;

R1a = 1./T1a;
R2a = 1./T2a;
R1b = R1a;

fb = 0.139 ; 
kba = 23;
kab = fb*kba;

M0b = fb;
M0a = 1;

% acquisition parameter
tp = 0.7*1e-3; % pulse
td = 2e-3;
Np = 100;

Gz = 15*1e-3; % gradient
r = 0;

% Different w1 values to compare [Hz]
w1_Hz_vals = [50, 100, 200, 500, 1000];
n_w1 = numel(w1_Hz_vals);

% dw sweep range
dw_Hz_vals = 1:10:5000;
n_dw = numel(dw_Hz_vals);

% Colors for different w1 (使用更鲜艳的配色)
colors = [0.0000, 0.4470, 0.7410;   % 蓝色
          0.8500, 0.3250, 0.0980;   % 橙色
          0.4660, 0.6740, 0.1880;   % 绿色
          0.4940, 0.1840, 0.5560;   % 紫色
          0.9290, 0.6940, 0.1250];  % 金色

%% Pre-compute common matrices
A_free = [-R2a,     0,        0,       0;
             0,  -R2a,        0,       0;
             0,     0,  -R1a-kab,      kba;
             0,     0,       kab,  -R1b-kba];

C = [0, 0, R1a*M0a, R1b*M0b]';

pa = 0 + gamma*Gz*r*td;

Rz = [cos(pa) , sin(pa), 0, 0;
      -sin(pa), cos(pa), 0, 0;
      0       , 0      , 1, 0;
      0       , 0      , 0, 1];

Mini = [0; 0; M0a; M0b];

%% Storage for results
M_BM = zeros(n_w1, n_dw);
M_RotMat = zeros(n_w1, n_dw);

%% Simulation for each w1
for iw = 1:n_w1
    w1 = 2*pi*w1_Hz_vals(iw);
    
    fprintf('Simulating w1 = %d Hz...\n', w1_Hz_vals(iw));
    
    for i = 1:n_dw
        dw = 2*pi*dw_Hz_vals(i);
        
        %% 1) BM solution
        M = Mini;
        Rrfb = RF_MT(T2b, w1, dw, 'SuperLorentzian');
        
        A = [-R2a,        dw,         0,       0;
              -dw,      -R2a,        w1,       0;
                0,       -w1,  -R1a-kab,      kba;
                0,         0,       kab,     -R1b-Rrfb-kba];
        
        for n = 1:Np
            M = expm(A.*tp)*(M + A\C) - A\C;
            M = expm(A_free.*td)*(M + A_free\C) - A_free\C;
            M = Rz*M;
        end
        M_BM(iw, i) = M(3);
        
        %% 2) Rotate Matrix Approximation
        M = Mini;
        
        E2a = exp(-(td+tp)/T2a);
        E1a = exp(-(td+tp)/T1a);
        E1b = exp(-R1b*(td+tp));
        expfac = exp(-(1+fb)*kba*(td+tp));
        
        E = diag([E2a E2a E1a E1b]);
        
        A_mat = [1+fb, 0, 0, 0;
                 0, 1+fb, 0, 0;
                 0, 0, 1+fb*expfac, 1-expfac;
                 0, 0, fb-fb*expfac, fb+expfac]/(1+fb);
        
        M0_vec = [0; 0; M0a*(1-E1a); M0b*(1-E1b)];
        
        Omega_eff = sqrt(w1^2 + dw^2);
        theta = Omega_eff * tp;
        alpha = w1 / Omega_eff;
        beta = dw / Omega_eff;
        
        R = [cos(theta),     -beta*sin(theta),  alpha*sin(theta), 0;
             beta*sin(theta), cos(theta)+(1-cos(theta))*alpha^2,   (1-cos(theta))*alpha*beta, 0;
            -alpha*sin(theta), (1-cos(theta))*alpha*beta, cos(theta)+(1-cos(theta))*beta^2, 0;
            0, 0, 0, exp(-Rrfb*tp)];
        
        for n = 1:Np
            M = Rz*R*A_mat*(E*M + M0_vec);
        end
        M_RotMat(iw, i) = M(3);
    end
end

%% Figure: BM彩色实线，RotMat深灰色虚线覆盖
fig = figure('Color', 'w', 'Position', [752, 451, 407, 344]);
hold on;

% 添加网格（先画网格，这样线条在上面）
grid on;
set(gca, 'GridLineStyle', '-', 'GridAlpha', 0.3, 'GridColor', [0.5 0.5 0.5]);

% 先画所有BM（彩色实线）
h_BM = gobjects(n_w1, 1);
for iw = 1:n_w1
    h_BM(iw) = plot(dw_Hz_vals, M_BM(iw, :), '-', 'LineWidth', 1.6, 'Color', colors(iw, :));
end

% 再画所有RotMat（深灰色虚线，覆盖在上面）
gray_color = [0.25, 0.25, 0.25];
for iw = 1:n_w1
    plot(dw_Hz_vals, M_RotMat(iw, :), '--', 'LineWidth', 1.6, 'Color', gray_color);
end

% 添加一条用于legend的灰色虚线
h_RotMat = plot(NaN, NaN, '--', 'LineWidth', 1.6, 'Color', gray_color);

% Legend
legend_str = cell(n_w1 + 1, 1);
for iw = 1:n_w1
    legend_str{iw} = sprintf('\\omega_1 = %d Hz', w1_Hz_vals(iw));
end
legend_str{n_w1 + 1} = 'Approx.';

lg = legend([h_BM; h_RotMat], legend_str, 'Location', 'southeast', 'FontSize', 10);
lg.FontWeight = font_weight;
lg.Box = 'off';  % 去掉legend边框，更简洁

% 坐标轴标签
xlabel('\Delta\omega (Hz)', 'FontSize', 13, 'FontWeight', font_weight);
ylabel('M_z', 'FontSize', 13, 'FontWeight', font_weight);

% 坐标轴范围
ylim([0, 1]);
xlim([0, 5000]);

% 设置坐标轴属性
ax = gca;
% ax.FontSize = 11;
% ax.FontWeight = font_weight;
% ax.LineWidth = 1;           % 坐标轴线条加粗
ax.Box = 'on';                % 显示完整边框
% ax.XMinorTick = 'on';         % 显示x轴小刻度
% ax.YMinorTick = 'on';         % 显示y轴小刻度
% ax.TickDir = 'out';           % 刻度向外
% ax.TickLength = [0.015, 0.015];

% 设置刻度值
ax.XTick = 0:500:5000;
ax.YTick = 0:0.2:1;

%% Print summary
fprintf('\n========== Summary ==========\n');
fprintf('tp = %.1f ms, td = %.1f ms, Np = %d\n', tp*1e3, td*1e3, Np);
fprintf('\nMax |BM - RotMat|:\n');
for iw = 1:n_w1
    max_diff = max(abs(M_BM(iw, :) - M_RotMat(iw, :)));
    fprintf('  w1 = %4d Hz: %.4f\n', w1_Hz_vals(iw), max_diff);
end

exportgraphics(gcf, 'figure matrix validation.png', 'Resolution', 600);