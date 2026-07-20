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

%% ==================== USER SETTINGS ====================
USE_BOLD_TEXT = true;

line_colors = [0.49 0.18 0.56;   % purple
            0.00 0.45 0.74;   % blue
            [255, 182, 193] / 255;
            0.47 0.67 0.19;   % green
            0.93 0.69 0.13;   % orange
            0.85 0.10 0.10;   % red
            [255, 247, 0] / 255; 
             [0.6902, 0.7686, 0.8706]];  

line_width       = 1.4;
title_font_size  = 14;
label_font_size  = 14;
axis_font_size   = 12;
legend_font_size = 12;

% --- Experimental operating point (vertical reference lines) ---
EXP_Gz_mTm = 25;    % experimental Gz [mT/m]  -> subplot (B)
EXP_Np     = 100;   % experimental Np         -> subplot (C)
exp_line_color = [0.15 0.15 0.15];
exp_line_style = ':';
exp_line_width = 1.6;
%% =======================================================

if USE_BOLD_TEXT, fontWeight = 'bold'; else, fontWeight = 'normal'; end

%% ---------------------- DANTE prep settings -----------------------------
gamma      = 2.6752e8;
tp         = 0.7e-3;
td         = 2e-3;
tg         = 1e-3;
phi        = 0;
dw         = 2*pi*1000;
w1_base    = 2*pi*100;
Gz_nominal = 25e-3;

v_arr_cm = [0, 0.1, 0.5, 1.0, 5.0, 10.0, 40.0, 80.0];
v_arr    = v_arr_cm * 1e-2;
n_v      = numel(v_arr);

angle_deg = linspace(0, 180, 181);
angle_rad = deg2rad(angle_deg);

Gz_mTm = linspace(0, 50, 101);
Gz_arr = Gz_mTm * 1e-3;

Np_max   = 1000;
Np_axis  = 0:Np_max;

Np_angle      = 100;
Np_grad       = 100;
angle_fix_rad = 0;
v_static      = 0;

%% ---------------------- TSE readout settings ----------------------------
TSE.N_startup       = 3;
TSE.ETL_acq         = 72;
TSE.ETL_tot         = TSE.N_startup + TSE.ETL_acq;
TSE.ESP             = 3.4e-3;
TSE.tp_ex           = 0.95e-3;
TSE.tp_ref          = 0.78e-3;
TSE.flip_ex         = pi/2;
TSE.flip_ref        = pi/2;
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
TSE.Gz_readout      = Gz_nominal;
TSE.tg_readout      = tg;

fprintf('TSE readout: ETL=%d (+%d startup), ESP=%.2f ms, TE_eff=%.2f ms\n', ...
        TSE.ETL_acq, TSE.N_startup, TSE.ESP*1e3, TSE.TEeff_ms);
fprintf('Flow during readout: %s   (Gz=%.1f mT/m, tg=%.2f ms)\n', ...
        string(TSE.flow_enable), TSE.Gz_readout*1e3, TSE.tg_readout*1e3);

%% ---------------------- Tissue definitions ------------------------------
WM    = struct('name','WM',   'T1',1084e-3,'T2',69e-3, 'T2b',10e-6, 'fb',0.139,'kba',23);
GM    = struct('name','GM',   'T1',1820e-3,'T2',99e-3, 'T2b',9.1e-6,'fb',0.05, 'kba',40);
Blood = struct('name','Blood','T1',1932e-3,'T2',275e-3,'T2b',280e-6,'fb',0.028,'kba',35);

%% ---------------------- WM static reference -----------------------------
[WMpar, WMrf] = build_tissue(WM, w1_base, dw, tp, TSE);

S_WM_angle = run_angle_sweep(Np_angle, angle_rad, v_static, ...
                             tp, td, tg, WMpar, WMrf, phi, Gz_nominal, TSE);
S_WM_grad  = run_gradient_sweep(Np_grad, Gz_arr, angle_fix_rad, v_static, ...
                                tp, td, tg, WMpar, WMrf, phi, TSE);
S_WM_Np    = run_Np_evolution(Np_max, angle_fix_rad, Gz_nominal, v_static, ...
                              tp, td, tg, WMpar, WMrf, phi, TSE);

%% ---------------------- GM static reference -----------------------------
[GMpar, GMrf] = build_tissue(GM, w1_base, dw, tp, TSE);

S_GM_angle = run_angle_sweep(Np_angle, angle_rad, v_static, ...
                             tp, td, tg, GMpar, GMrf, phi, Gz_nominal, TSE);
S_GM_grad  = run_gradient_sweep(Np_grad, Gz_arr, angle_fix_rad, v_static, ...
                                tp, td, tg, GMpar, GMrf, phi, TSE);
S_GM_Np    = run_Np_evolution(Np_max, angle_fix_rad, Gz_nominal, v_static, ...
                              tp, td, tg, GMpar, GMrf, phi, TSE);

%% ---------------------- Blood ------------------------------------------
[Bpar, Brf] = build_tissue(Blood, w1_base, dw, tp, TSE);

%% ---------------------- Figure ------------------------------------------
bigFig = figure('Color','w', ...
                'Name','DANTE-TSE (flow-aware) - Blood vs Static WM/GM', ...
                'Position',[652         410        1159         400]);
tl = tiledlayout(bigFig,1,3,'TileSpacing','compact','Padding','compact');
tile_labels = {'(A)','(B)','(C)'};

%% --- (A) Flow direction sweep ---
ax1 = nexttile(tl,1); hold(ax1,'on'); grid(ax1,'on');

plot(ax1, angle_deg, S_WM_angle, '--','LineWidth',line_width, ...
     'Color',[0.5 0.5 0.5],'DisplayName','Static WM');
plot(ax1, angle_deg, S_GM_angle, '--','LineWidth',line_width, ...
     'Color',[0.2 0.6 0.2],'DisplayName','Static GM');

for j = 1:n_v
    Sig_a = run_angle_sweep(Np_angle, angle_rad, v_arr(j), ...
                            tp, td, tg, Bpar, Brf, phi, Gz_nominal, TSE);
    plot(ax1, angle_deg, Sig_a, 'LineWidth',line_width,'Color',line_colors(j,:));
end

xlabel(ax1,'\theta (degrees)','FontWeight',fontWeight);
ax1.XLabel.FontSize = label_font_size;
ylabel(ax1,'M_{xy}','FontWeight',fontWeight);
ax1.YLabel.FontSize = label_font_size;
title(ax1,sprintf('%s Flow Direction',tile_labels{1}),'FontWeight',fontWeight);
ax1.Title.FontSize = title_font_size;
ax1.XLim = [0,180]; ax1.XTick = 0:30:180;
ax1.YLim(1) = -0.2;
ax1.FontSize = axis_font_size;

%% --- (B) Gradient sweep ---
ax2 = nexttile(tl,2); hold(ax2,'on'); grid(ax2,'on');

plot(ax2, Gz_mTm, S_WM_grad, '--','LineWidth',line_width, ...
     'Color',[0.5 0.5 0.5],'DisplayName','Static WM');
plot(ax2, Gz_mTm, S_GM_grad, '--','LineWidth',line_width, ...
     'Color',[0.2 0.6 0.2],'DisplayName','Static GM');

for j = 1:n_v
    Sig_g = run_gradient_sweep(Np_grad, Gz_arr, angle_fix_rad, v_arr(j), ...
                               tp, td, tg, Bpar, Brf, phi, TSE);
    plot(ax2, Gz_mTm, Sig_g, 'LineWidth',line_width,'Color',line_colors(j,:));
end

% --- Vertical reference line @ experimental Gz ---
xline(ax2, EXP_Gz_mTm, exp_line_style, ...
      sprintf('G_z = %g mT/m', EXP_Gz_mTm), ...
      'Color', exp_line_color, 'LineWidth', exp_line_width, ...
      'LabelOrientation','horizontal','LabelVerticalAlignment','top', ...
      'LabelHorizontalAlignment','center', ...
      'FontWeight',fontWeight,'FontSize',axis_font_size, ...
      'HandleVisibility','off');

xlabel(ax2,'G_z (mT/m)','FontWeight',fontWeight);
ax2.XLabel.FontSize = label_font_size;
ylabel(ax2,'');
title(ax2,sprintf('%s Gradient',tile_labels{2}),'FontWeight',fontWeight);
ax2.Title.FontSize = title_font_size;
ax2.XLim = [min(Gz_mTm), max(Gz_mTm)];
ax2.FontSize = axis_font_size;

%% --- (C) Np evolution ---
ax3 = nexttile(tl,3); hold(ax3,'on'); grid(ax3,'on');

plot(ax3, Np_axis, S_WM_Np, '--','LineWidth',line_width, ...
     'Color',[0.5 0.5 0.5],'DisplayName','Static WM');
plot(ax3, Np_axis, S_GM_Np, '--','LineWidth',line_width, ...
     'Color',[0.2 0.6 0.2],'DisplayName','Static GM');

for j = 1:n_v
    Sig_n = run_Np_evolution(Np_max, angle_fix_rad, Gz_nominal, v_arr(j), ...
                             tp, td, tg, Bpar, Brf, phi, TSE);
    plot(ax3, Np_axis, Sig_n, 'LineWidth',line_width,'Color',line_colors(j,:));
end

% --- Vertical reference line @ experimental Np ---
xline(ax3, EXP_Np, exp_line_style, ...
      sprintf('N_p = %g', EXP_Np), ...
      'Color', exp_line_color, 'LineWidth', exp_line_width, ...
      'LabelOrientation','horizontal','LabelVerticalAlignment','top', ...
      'LabelHorizontalAlignment','center', ...
      'FontWeight',fontWeight,'FontSize',axis_font_size, ...
      'HandleVisibility','off');

xlabel(ax3,'N_p','FontWeight',fontWeight);
ax3.XLabel.FontSize = label_font_size;
ylabel(ax3,'');
title(ax3,sprintf('%s Number of Pulse',tile_labels{3}),'FontWeight',fontWeight);
ax3.Title.FontSize = title_font_size;
ax3.XLim = [0, Np_max];
ax3.FontSize = axis_font_size;

legend_entries = [{'Static WM','Static GM'}, ...
                  arrayfun(@(vc) sprintf('v = %.1f cm/s',vc), v_arr_cm, ...
                           'UniformOutput',false)];
leg = legend(ax3, legend_entries, 'Location','east','Box','off');
leg.FontSize = legend_font_size;

% sgtitle(sprintf(['DANTE + TSE readout (flow-aware) |  TE_{eff} = %.1f ms  ' ...
%                 '(echo #%d, ETL = %d+%d, ESP = %.1f ms, refoc = 90°)'], ...
%                 TSE.TEeff_ms, TSE.echo_idx_TEeff, TSE.N_startup, ...
%                 TSE.ETL_acq, TSE.ESP*1e3), ...
%         'FontWeight',fontWeight,'FontSize',12);

exportgraphics(gcf,'figure parameter selection1_revision.png','Resolution',600);
fprintf('Done.\n');


%% ========================================================================
%                               Local functions
% ========================================================================
function [par, rf] = build_tissue(T, w1_base, dw, tp, TSE)
    par.T1=T.T1; par.T2=T.T2; par.T2b=T.T2b;
    par.fb=T.fb; par.kba=T.kba; par.kab=T.fb*T.kba;
    par.R1a=1/T.T1; par.R2a=1/T.T2; par.R1b=par.R1a;
    par.M0a=1; par.M0b=T.fb;

    Omega_eff = sqrt(w1_base^2 + dw^2);
    rf.beta_dante  = Omega_eff*tp;
    rf.theta_dante = acos(max(min(dw/Omega_eff,1),-1));
    rf.Rrfb_dante  = RF_MT(T.T2b, w1_base, dw, 'SuperLorentzian');

    rf.Rrfb_ex  = RF_MT(T.T2b, TSE.w1_ex,  TSE.dw_tse, 'SuperLorentzian');
    rf.Rrfb_ref = RF_MT(T.T2b, TSE.w1_ref, TSE.dw_tse, 'SuperLorentzian');
end

% ------------------------------------------------------------------------
function S = run_DANTE_block(Np, v, angle, tp, td, tg, par, rf, phi, Gz)
    S = [0; 0; par.M0a; par.M0b];
    for n = 1:Np
        S = epgmt_RF(S, rf.theta_dante, rf.beta_dante, phi, rf.Rrfb_dante, tp);
        S = epgmt_Ex_Relax(S, par.R1a, par.R2a, par.R1b, par.kab, par.kba, ...
                           tp, td, par.M0a, par.M0b);
        S = epgmt_Grad(S);
        S = epgmt_Flow(S, Gz, tg, v, angle);
    end
end

% ------------------------------------------------------------------------
function sig = run_TSE_readout(S, par, rf, TSE, v, flow_angle)
    if nargin < 5, v = 0;            end
    if nargin < 6, flow_angle = 0;   end
    do_flow = TSE.flow_enable && (v ~= 0);

    S = epgmt_RF(S, TSE.theta_tse, TSE.flip_ex, TSE.phi_ex, rf.Rrfb_ex, TSE.tp_ex);
    S = epgmt_Ex_Relax(S, par.R1a, par.R2a, par.R1b, par.kab, par.kba, ...
                       0, TSE.gap_ex_to_ref1, par.M0a, par.M0b);
    if do_flow
        S = epgmt_Flow(S, TSE.Gz_readout, TSE.tg_readout, v, flow_angle);
    end

    sig = NaN;
    for e = 1:TSE.ETL_tot
        if e > 1
            S = epgmt_Ex_Relax(S, par.R1a, par.R2a, par.R1b, par.kab, par.kba, ...
                               0, TSE.half_esp, par.M0a, par.M0b);
            if do_flow
                S = epgmt_Flow(S, TSE.Gz_readout, TSE.tg_readout, v, flow_angle);
            end
        end
        S = epgmt_Grad(S);
        S = epgmt_RF(S, TSE.theta_tse, TSE.flip_ref, TSE.phi_ref, ...
                     rf.Rrfb_ref, TSE.tp_ref);
        S = epgmt_Grad(S);
        S = epgmt_Ex_Relax(S, par.R1a, par.R2a, par.R1b, par.kab, par.kba, ...
                           0, TSE.half_esp, par.M0a, par.M0b);
        if do_flow
            S = epgmt_Flow(S, TSE.Gz_readout, TSE.tg_readout, v, flow_angle);
        end
        if e == TSE.echo_idx_TEeff
            sig = real(S(1,1));
            return;
        end
    end
end

% ------------------------------------------------------------------------
function Sig = run_angle_sweep(Np, angle_rad, v, tp, td, tg, par, rf, phi, Gz, TSE)
    Sig = zeros(1, numel(angle_rad));
    for k = 1:numel(angle_rad)
        ang  = angle_rad(k);
        S    = run_DANTE_block(Np, v, ang, tp, td, tg, par, rf, phi, Gz);
        Sig(k) = run_TSE_readout(S, par, rf, TSE, v, ang);
    end
end

% ------------------------------------------------------------------------
function Sig = run_gradient_sweep(Np, Gz_arr, angle_fix, v, tp, td, tg, par, rf, phi, TSE)
    Sig = zeros(1, numel(Gz_arr));
    for k = 1:numel(Gz_arr)
        S = run_DANTE_block(Np, v, angle_fix, tp, td, tg, par, rf, phi, Gz_arr(k));
        Sig(k) = run_TSE_readout(S, par, rf, TSE, v, angle_fix);
    end
end

% ------------------------------------------------------------------------
function Sig_evo = run_Np_evolution(Np_max, angle_fix, Gz, v, tp, td, tg, par, rf, phi, TSE)
    Sig_evo = zeros(1, Np_max+1);

    S0 = [0; 0; par.M0a; par.M0b];
    Sig_evo(1) = run_TSE_readout(S0, par, rf, TSE, v, angle_fix);

    S = S0;
    for n = 1:Np_max
        S = epgmt_RF(S, rf.theta_dante, rf.beta_dante, phi, rf.Rrfb_dante, tp);
        S = epgmt_Ex_Relax(S, par.R1a, par.R2a, par.R1b, par.kab, par.kba, ...
                           tp, td, par.M0a, par.M0b);
        S = epgmt_Grad(S);
        S = epgmt_Flow(S, Gz, tg, v, angle_fix);

        Sig_evo(n+1) = run_TSE_readout(S, par, rf, TSE, v, angle_fix);
    end
end