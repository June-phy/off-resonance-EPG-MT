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

% 自定义线条颜色
% 方法1: 使用颜色字符
% line_colors = {'b', 'r', 'g', 'm', 'c', 'k'};  % 蓝、红、绿、品红、青、黑

% 方法2: 使用RGB值 (取消注释以使用)
% line_colors = [0.0000, 0.4470, 0.7410;   % 蓝色
%                0.8500, 0.3250, 0.0980;   % 橙色
%                0.9290, 0.6940, 0.1250;   % 黄色
%                0.4940, 0.1840, 0.5560;   % 紫色
%                0.4660, 0.6740, 0.1880;   % 绿色
%                0.3010, 0.7450, 0.9330];  % 浅蓝色

line_colors = [0.0000, 0.4470, 0.7410;   % 蓝色
               0.8500, 0.3250, 0.0980;   % 橙色
               0.9290, 0.6940, 0.1250;   % 黄色
               0.4940, 0.1840, 0.5560;   % 紫色
               0.4660, 0.6740, 0.1880;   % 绿色
               0.3010, 0.7450, 0.9330];  % 浅蓝色

% 线宽设置
line_width = 1.5;
%% =======================================================

%% Constants & shared settings
gamma  = 2.6752e8;      % rad/(s*T)
tp     = 0.7e-3;        % [s] RF pulse duration
td     = 2.0e-3;        % [s] spacing between DANTE modules
Np     = 350;           % number of pulses
pulse_index = 0:Np;     % x-axis vector for plots

Gz     = 25e-3;          % [T/m] gradient amplitude
tg     = 1e-3;          % [s] gradient duration (used in flow phase term)

w1     = 2*pi*100;      % [rad/s] RF amplitude
dw     = 2*pi*1000;     % [rad/s] off-resonance
phi    = 0;             % [rad] RF phase (can be extended to phase cycling)
angle  = 0;             % [rad] flow direction (0 => aligned with gradient)

v_arr  = [0, 0.1, 0.5, 1.0, 5.0, 10]*1e-2;   % [m/s] speeds (converted from cm/s)
n_v    = numel(v_arr);

N_iso  = 1000;          % number of isochromats in distributed simulation
L      = 0.01;          % [m] spatial extent for isochromat distribution
pos_arr = linspace(-L/2, L/2, N_iso);

%% Parameter sets for white matter & blood
% White matter parameters
paramSets(1).label = 'White Matter';
paramSets(1).T1    = 1084e-3;
paramSets(1).T2    = 69e-3;
paramSets(1).T2b   = 10e-6;
paramSets(1).fb    = 0.139;
paramSets(1).kba   = 23;

% Blood parameters
paramSets(2).label = 'Blood';
paramSets(2).T1    = 1932e-3;
paramSets(2).T2    = 275e-3;
paramSets(2).T2b   = 280e-6;
paramSets(2).fb    = 0.028;
paramSets(2).kba   = 35;

%% Containers for results
for p = 1:numel(paramSets)
    results(p).Mz_single_iso  = zeros(n_v, Np+1);
    results(p).Mz_multi_iso   = zeros(n_v, Np+1);
    results(p).Mz_epg         = zeros(n_v, Np+1);
end

%% Main simulation loops
for p = 1:numel(paramSets)
    % --- Parameter unpack ---
    T1   = paramSets(p).T1;
    T2   = paramSets(p).T2;
    T2b  = paramSets(p).T2b;
    fb   = paramSets(p).fb;
    kba  = paramSets(p).kba;
    kab  = fb * kba;

    R1a = 1/T1;
    R2a = 1/T2;
    R1b = R1a;

    M0a = 1;
    M0b = fb;
    Mini = [0; 0; M0a; M0b];

    % Effective RF properties
    Rrfb = RF_MT(T2b, w1, dw, 'SuperLorentzian');
    Omega_eff = sqrt(w1^2 + dw^2);
    beta  = Omega_eff * tp;
    theta = acos(dw / Omega_eff);

    % Bloch-McConnell matrices
    A = [ -R2a,      dw,       0,               0;
          -dw,      -R2a,      w1,              0;
           0,       -w1,    -R1a-kab,          kba;
           0,        0,        kab,    -R1b-Rrfb-kba ];

    A_free = [ -R2a,    0,        0,         0;
                0,    -R2a,       0,         0;
                0,      0,    -R1a-kab,     kba;
                0,      0,       kab,    -R1b-kba ];

    C = [0; 0; R1a*M0a; R1b*M0b];

    % Precompute propagators & steady-state offsets
    pulseProp     = expm(A * tp);
    freeProp      = expm(A_free * td);
    Ainv_C        = A \ C;
    Afree_inv_C   = A_free \ C;

    % --- Loop over velocities ---
    for j = 1:n_v
        v = v_arr(j);

        %% (1) Single isochromat (Bloch-McConnell integration with flow phase)
        
        M = Mini;
        Mz_single = zeros(1, Np+1);
        Mz_single(1) = Mini(3);

        for n = 1:Np
            M = pulseProp * (M + Ainv_C) - Ainv_C;
            M = freeProp  * (M + Afree_inv_C) - Afree_inv_C;

            pa = n * gamma * Gz * v * tg^2; % accumulated phase
            rotZ = [ cos(pa),  sin(pa), 0, 0;
                    -sin(pa),  cos(pa), 0, 0;
                         0,        0,  1, 0;
                         0,        0,  0, 1 ];

            M = rotZ * M;
            Mz_single(n+1) = real(M(3));
        end
        results(p).Mz_single_iso(j, :) = Mz_single;

        %% (2) Distributed isochromats (1000 positions averaged)
        tic;
        M_all = repmat(Mini, 1, N_iso);
        Mz_total = zeros(1, Np+1);
        Mz_total(1) = Mini(3);

        for n = 1:Np
            for iso = 1:N_iso
                pos = pos_arr(iso);

                Mtmp = pulseProp * (M_all(:, iso) + Ainv_C) - Ainv_C;
                Mtmp = freeProp  * (Mtmp + Afree_inv_C)   - Afree_inv_C;

                pa = gamma * Gz * pos * tg ...
                   - gamma * Gz * v * tg^2 / 2 ...
                   + n * gamma * Gz * v * tg^2;

                rotZ = [ cos(pa),  sin(pa), 0, 0;
                        -sin(pa),  cos(pa), 0, 0;
                             0,        0,  1, 0;
                             0,        0,  0, 1 ];
                M_all(:, iso) = rotZ * Mtmp;
            end
            Mz_total(n+1) = mean(real(M_all(3, :)));
        end
        results(p).Mz_multi_iso(j, :) = Mz_total;
        time_1000_iso = toc; % Stop timing for 1000 isochromats
        fprintf('1000 Isochromats (%s, v = %.2f cm/s): %.6f s\n', ...
                paramSets(p).label, v*100, time_1000_iso);

        %% (3) EPG simulation (discrete states)
        tic;
        S = [0; 0; M0a; M0b];
        Mz_epg = zeros(1, Np+1);
        Mz_epg(1) = M0a;

        for n = 1:Np
            S = epgmt_RF(S, theta, beta, phi, Rrfb, tp);
            S = epgmt_Ex_Relax(S, R1a, R2a, R1b, kab, kba, tp, td, M0a, M0b);
            S = epgmt_Grad(S);
            S = epgmt_Flow(S, Gz, tg, v, angle);

            Mz_epg(n+1) = real(S(3, 1));
        end
        results(p).Mz_epg(j, :) = Mz_epg;
        time_epg = toc; % Stop timing for EPG simulation
        fprintf('EPG Simulation (%s, v = %.2f cm/s): %.6f s\n', ...
                paramSets(p).label, v*100, time_epg);
    end
end

%% Plot montage (2 rows x 3 columns)
velocity_labels = {'v = 0.0 cm/s', 'v = 0.1 cm/s', 'v = 0.5 cm/s', ...
                   'v = 1.0 cm/s', 'v = 5.0 cm/s', 'v = 10.0 cm/s'};

subplot_labels = {'(A)', '(B)', '(C)', '(D)', '(E)', '(F)'};

% Set FontWeight based on toggle
if USE_BOLD_TEXT
    fontWeight = 'bold';
else
    fontWeight = 'normal';
end

figure('Color', 'w', 'Position', [100, 100, 1200, 700]);
t = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

for p = 1:numel(paramSets)
    % --- Single isochromat ---
    tile_idx = (p-1)*3 + 1;
    ax = nexttile(tile_idx); hold(ax, 'on'); grid(ax, 'on');
    for j = 1:n_v
        % 根据 line_colors 的类型选择颜色
        if iscell(line_colors)
            plot(ax, pulse_index, results(p).Mz_single_iso(j, :), ...
                 'LineWidth', line_width, 'Color', line_colors{j});
        else
            plot(ax, pulse_index, results(p).Mz_single_iso(j, :), ...
                 'LineWidth', line_width, 'Color', line_colors(j, :));
        end
    end
    title(ax, sprintf('%s %s – Single BM Isochromat', subplot_labels{tile_idx}, paramSets(p).label), ...
          'FontWeight', fontWeight);
    ylabel(ax, 'M_z', 'FontWeight', fontWeight);
    if p == numel(paramSets)
        xlabel(ax, 'N_p', 'FontWeight', fontWeight);
    end
    if p == 1
        legend(ax, velocity_labels, 'Location', 'southwest', 'Box', 'off');
    end

    % --- 1000 isochromats ---
    tile_idx = (p-1)*3 + 2;
    ax = nexttile(tile_idx); hold(ax, 'on'); grid(ax, 'on');
    for j = 1:n_v
        if iscell(line_colors)
            plot(ax, pulse_index, results(p).Mz_multi_iso(j, :), ...
                 'LineWidth', line_width, 'Color', line_colors{j});
        else
            plot(ax, pulse_index, results(p).Mz_multi_iso(j, :), ...
                 'LineWidth', line_width, 'Color', line_colors(j, :));
        end
    end
    title(ax, sprintf('%s %s – 1000 BM Isochromats', subplot_labels{tile_idx}, paramSets(p).label), ...
          'FontWeight', fontWeight);
    ylabel(ax, '', 'FontWeight', fontWeight);
    if p == numel(paramSets)
        xlabel(ax, 'N_p', 'FontWeight', fontWeight);
    end

    % --- EPG ---
    tile_idx = (p-1)*3 + 3;
    ax = nexttile(tile_idx); hold(ax, 'on'); grid(ax, 'on');
    for j = 1:n_v
        if iscell(line_colors)
            plot(ax, pulse_index, results(p).Mz_epg(j, :), ...
                 'LineWidth', line_width, 'Color', line_colors{j});
        else
            plot(ax, pulse_index, results(p).Mz_epg(j, :), ...
                 'LineWidth', line_width, 'Color', line_colors(j, :));
        end
    end
    title(ax, sprintf('%s Off resonance EPG-MT', subplot_labels{tile_idx}), ...
          'FontWeight', fontWeight);
    ylabel(ax, '', 'FontWeight', fontWeight);
    if p == numel(paramSets)
        xlabel(ax, 'N_p', 'FontWeight', fontWeight);
    end
end

exportgraphics(gcf, 'figure epg validation.png', 'Resolution', 600);