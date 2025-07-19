%% ========================================================================
%  運動方程式の項の比較 (慣性項 vs JxB力) - トロイダル成分 2Dプロット
% =========================================================================
clear;
close all;
clc;

%% ★★★ ユーザー設定 ★★★
% --- プロットするケースを選択 ---
% Case-I
date = 240610;
% JxB, gridデータ
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240610023.mat');
% ne, Teデータ
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240610\triple_data2D_case-I.mat');
% vt (Flow) データ
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\MachProbe_data\mat\240610\18-47_Iratio.mat'); % Flow計算済みのMATファイルを指定

% % case-O
% date = 240611;
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240611055.mat');
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240611\triple_data2D_case-O.mat');
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\MachProbe_data\mat\240611\47-70_Iratio.mat'); % Flow計算済みのMATファイルを指定

% --- プロットのタイミングを設定 ---
start_time_us = 480;
time_step_us = 5; % 時間間隔を少し広めに設定
num_plots_y = 2;
num_plots_x = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- 1. データ準備と慣性項 (dv_t/dt) の事前計算 ---
fprintf('Calculating toroidal acceleration...\n');
% --- 各データの時間・空間軸を取得
jxb_time_axis = data2D.trange;
ne_time_axis = time_values;
ne_R_axis = R_values;
vt_time_axis = MPdata2D.trange;
vt_R_axis = MPdata2D.rq(:,1);
vt_Z_axis = MPdata2D.zq(1,:);

% --- トロイダル速度 vt [m/s] を取得 ---
% MP_calcで計算されたflow_forplotがvtに相当
% 次元: (time, r, z)
vt_3D = squeeze(MPdata2D.flow_forplot);

% --- トロイダル加速度 at = dv_t/dt [m/s^2] を計算 ---
% gradient関数で時間微分
time_step_s = (vt_time_axis(2) - vt_time_axis(1)) * 1e-6;
at_3D = gradient(vt_3D, time_step_s); % 時間軸(1次元目)に沿って微分

mi = 1.672e-27; % イオン質量 [kg]
ne_data = squeeze(triple_data2D.ne_avg); % (R, Time)

%% --- 2. プロット処理 ---
figure('WindowState', 'maximized');
total_plots = num_plots_y * num_plots_x;

for plot_idx = 1:total_plots
    time_in_us = start_time_us + (plot_idx - 1) * time_step_us;
    
    subplot(num_plots_y, num_plots_x, plot_idx);
    
    % --- a) JxB力のトロイダル成分を計算 ---
    [~, target_time_idx_jxb] = min(abs(jxb_time_axis - time_in_us));
    Jr = data2D.Jr(:, :, target_time_idx_jxb); Jt = data2D.Jt(:, :, target_time_idx_jxb); Jz = data2D.Jz(:, :, target_time_idx_jxb);
    Br = data2D.Br(:, :, target_time_idx_jxb); Bt = data2D.Bt(:, :, target_time_idx_jxb); Bz = data2D.Bz(:, :, target_time_idx_jxb);
    J_vectors = [Jr(:), Jt(:), Jz(:)];
    B_vectors = [Br(:), Bt(:), Bz(:)];
    j_B_force_vectors = cross(J_vectors, B_vectors, 2);
    % トロイダル成分(2番目)を2Dマップに戻す
    j_B_t_2D = reshape(j_B_force_vectors(:, 2), size(Jr));

    % --- b) 慣性項の計算 ---
    % 1. 現在時刻の加速度マップ at(r,z) を取得
    [~, target_time_idx_vt] = min(abs(vt_time_axis - time_in_us));
    at_2D = squeeze(at_3D(target_time_idx_vt, :, :)); % (r, z)
    
    % 2. 現在時刻の密度プロファイル ne(r) を取得
    [~, ne_time_idx] = min(abs(ne_time_axis - time_in_us));
    ne_profile_1D = ne_data(:, ne_time_idx);
    
    % 3. 密度プロファイルをvtのグリッドに合わせて2Dマップ化 (z方向に一様と仮定)
    ne_interp_1D = interp1(ne_R_axis, ne_profile_1D, vt_R_axis, 'pchip', 'extrap');
    ne_2D_aligned = repmat(ne_interp_1D, 1, length(vt_Z_axis));
    
    % 4. 慣性項の力 F_inertia = ni*mi*at を計算
    inertial_force_t_2D = (ne_2D_aligned * mi) .* at_2D;
    
    % --- c) プロット ---
    % JxB力をカラーコンターでプロット
    contourf(grid2D.zq, grid2D.rq, j_B_t_2D / 1000, 50, 'LineStyle', 'none');
    
    hold on;
    % 慣性項を線コンターで重ねてプロット
    contour(grid2D.zq, grid2D.rq, inertial_force_t_2D / 1000, 10, 'r--', 'LineWidth', 1);
    
    % 磁気面(psi=0)を太線でプロット
    contour(grid2D.zq, grid2D.rq, data2D.psi(:,:,target_time_idx_jxb), [0 0], 'k', 'LineWidth', 2);
    hold off;

    % --- グラフ装飾 ---
    h = colorbar;
    ylabel(h, '[kN/m^3]');
    clim([-100 100]); % カラーマップの範囲を調整
    xlabel('Z [m]');
    ylabel('R [m]');
    title(['t = ', num2str(time_in_us), ' us']);
    axis equal;
    
end
sgtitle(['Toroidal Force Analysis (J×B)_t vs m_in_ia_t for Date: ', num2str(date)], 'FontSize', 16, 'FontWeight', 'bold');