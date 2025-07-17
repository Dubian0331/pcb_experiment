clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% memo
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case-I:6/10, Xpoint:-0.01m -> val(;, 20, target_time) & val(:, 7, target_time), Target time: "400 + target_time". 499us
% ... (既存のmemoコメントはそのまま) ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ★★★ ユーザー設定 ★★★
% --- プロットするケースを選択 ---
% Case-I
date = '240610';
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240610018.mat');
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240610\grad_P_data_240610.mat');
z_index = 10; % プロットするzの位置 (z=0に対応)

% % case-O
% date = '240611';
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240611055.mat');
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\grad_P_data_240611.mat');
% z_index = 12; % プロットするzの位置 (z=0に対応)

% --- プロットのタイミングを設定 ---
start_time_us = 480; % プロットを開始する時刻 [us]
time_step_us = 1;    % プロットの時間間隔 [us]
num_plots_y = 5;     % 縦に並べるプロット数
num_plots_x = 6;     % 横に並べるプロット数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% --- プロット処理 ---
figure('WindowState', 'maximized');
% グリッドデータ (r方向)
r = grid2D.rq(:,1); % r方向の座標軸

% 指定されたタイミングでループ処理
total_plots = num_plots_y * num_plots_x;
for plot_idx = 1:total_plots
    % 現在のプロット時刻を計算
    time_in_us = start_time_us + (plot_idx - 1) * time_step_us;
    % data2Dのインデックスに変換
    target_time = time_in_us - 399;
    
    % サブプロットを作成
    subplot(num_plots_y, num_plots_x, plot_idx);
    hold on;
    grid on;
    
    % --- JxBの計算 (既存ロジック) ---
    Jt = data2D.Jt(:, :, target_time);
    Bz = data2D.Bz(:, :, target_time);
    j_B_r = Jt .* Bz; % r成分 (Jt*Bz - Jz*Bt のうち、JtBzのみ)
    jxb_to_plot = j_B_r(:, z_index) / 1000; % [kN/m^3] に変換
    % jxb_to_plot = j_B_r(:, z_index);
    
    plot(r, jxb_to_plot, 'b-', 'LineWidth', 1.5, 'DisplayName', 'J_tB_z');
    
    % --- 対応する時刻の圧力勾配を取得・補間 ---
    [~, grad_P_time_idx] = min(abs(time_axis - time_in_us));
    grad_P_to_plot = grad_P_2D(grad_P_time_idx, :);
    grad_P_aligned = interp1(R_axis, -grad_P_to_plot, r, 'pchip', 'extrap') / 1000; % [kN/m^3] に変換
    
    plot(r, grad_P_aligned, 'r--', 'LineWidth', 1.5, 'DisplayName', '-grad P');
    
    % --- ★★★ 追加: JxB と -grad P の和をプロット ★★★ ---
    force_sum = jxb_to_plot + grad_P_aligned;
    % plot(r, force_sum, 'k-', 'LineWidth', 2, 'DisplayName', 'Sum');
    
    % --- グラフ装飾 ---
    yline(0, 'k-', 'HandleVisibility', 'off');
    title(['t = ', num2str(time_in_us), ' us']);
    legend('FontSize', 7, 'Location', 'southwest');
    
    % 軸ラベル (外側のみに表示)
    if plot_idx > total_plots - num_plots_x
        xlabel('r [m]');
    end
    if mod(plot_idx - 1, num_plots_x) == 0
        ylabel('Force Density [kN/m^{3}]');
    end
    xlim([0.1, 0.3]);
end

% 全体にタイトルを追加
sgtitle(['Force Balance Analysis (J_tB_z vs -grad P) for Date: ', date], 'FontSize', 16, 'FontWeight', 'bold');

fprintf("save your file %s\n", date);