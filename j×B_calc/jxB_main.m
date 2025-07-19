%%%%%%%%%%%%%%%%%%%%%%%%%%%
% memo
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case-I:6/10, Xpoint:-0.01m -> val(;, 20, target_time) & val(:, 7, target_time), Target time: "400 + target_time". 499us
% startFile = 4289; 
% endFile = 4319;   
% X_point_is_here = 0.24;
% target_time = 99;

% Case-O:6/11, Xpoint:0.03m -> val(;, 22, target_time) & val(:, 9, target_time), file_start:4363, file_end:4389
% startFile = 4363; % 開始ファイル番号
% endFile = 4389;   % 終了ファイル番号
% X_point_is_here = 0.20;
% target_time = 90;
% target_time = 112;

%% ★★★ ユーザー設定 ★★★
% --- プロットするケースを選択 ---
% Case-I
% date = 240610;
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240610018.mat');
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240610\grad_P_data_240610.mat');
% z_index = 10; % プロットするzの位置 (z=0に対応)

% % case-O
date = 240611;
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240611055.mat');
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240611\grad_P_data_240611.mat');
z_index = 11; % プロットするzの位置 (z=0に対応)

% --- プロットのタイミングを設定 ---
start_time_us = 480; % プロットを開始する時刻 [us]
time_step_us = 1;    % プロットの時間間隔 [us]
num_plots_y = 5;     % 縦に並べるプロット数
num_plots_x = 6;     % 横に並べるプロット数
isSave = false; % 画像を保存するかどうか
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- プロット処理 ---
figure('WindowState', 'maximized');
r_axis_jxb = grid2D.rq(:,1); % JxBのR座標軸 (r)
jxb_time_axis = data2D.trange; % JxBの時間軸 (t)

total_plots = num_plots_y * num_plots_x;
for plot_idx = 1:total_plots
    time_in_us = start_time_us + (plot_idx - 1) * time_step_us;
    
    subplot(num_plots_y, num_plots_x, plot_idx);
    hold on;
    grid on;
    
    % --- JxBの全方向計算 ---
    % 1. 指定した時刻に最も近いインデックスを探す
    [~, target_time] = min(abs(jxb_time_axis - time_in_us));

    % 2. 必要なJとBの全成分データを取得 (次元は r, z)
    Jr = data2D.Jr(:, :, target_time);
    Jt = data2D.Jt(:, :, target_time);
    Jz = data2D.Jz(:, :, target_time);
    Br = data2D.Br(:, :, target_time);
    Bt = data2D.Bt(:, :, target_time);
    Bz = data2D.Bz(:, :, target_time);

    % 3. ベクトル化して外積を効率的に計算
    J_vectors = [Jr(:), Jt(:), Jz(:)];
    B_vectors = [Br(:), Bt(:), Bz(:)];
    j_B_force_vectors = cross(J_vectors, B_vectors, 2);

    % 4. r方向(半径方向)の成分を抜き出す
    j_B_r_raw = j_B_force_vectors(:, 1); % 1番目の列がr成分
    j_B_r_2D = reshape(j_B_r_raw, size(Jr));
    % jt×Bzのみ
    j_t_B_z_term = (Jt .* Bz) / 1000; % [kN/m^3]
    jtbz_term_to_plot = j_t_B_z_term(:, z_index);
    % jz×Btのみ
    j_z_B_t_term = (-Jz .* Bt) / 1000;
    jzb_term_to_plot = j_z_B_t_term(:, z_index);

    % 5. プロット用にzの位置でスライスし、単位を変換
    jxb_r_to_plot = j_B_r_2D(:, z_index) / 1000; % [kN/m^3] に変換

    % 6. r方向成分をプロット
    plot(r_axis_jxb, jxb_r_to_plot, 'b-', 'LineWidth', 1.5, 'DisplayName', '(J×B)_r');
    
    % --- grad Pのプロット---
    [~, grad_P_time_idx] = min(abs(time_axis - time_in_us));
    grad_P_to_plot = grad_P_2D(grad_P_time_idx, :);
    grad_P_aligned = interp1(R_axis, -grad_P_to_plot, r_axis_jxb, 'pchip', 'extrap') / 1000;
    
    % plot(r_axis_jxb, grad_P_aligned, 'r--', 'LineWidth', 1.5, 'DisplayName', '-∇P_r');

    plot(r_axis_jxb, jtbz_term_to_plot, 'g:', 'LineWidth', 2, 'DisplayName', 'J_tB_z term');

    plot(r_axis_jxb, jzb_term_to_plot, 'm-.', 'LineWidth', 2, 'DisplayName', '-J_zB_t term');
    
    % --- 和のプロット (力の釣り合いを確認) ---
    force_sum = jxb_r_to_plot + grad_P_aligned;
    % plot(r_axis_jxb, force_sum, 'k-', 'LineWidth', 2, 'DisplayName', 'Sum');
    
    % --- グラフ装飾 ---
    yline(0, 'k-', 'HandleVisibility', 'off');
    title(['t = ', num2str(time_in_us), ' us']);
    if plot_idx == 1, legend('FontSize', 7); end
    if plot_idx > total_plots - num_plots_x, xlabel('r [m]'); end
    if mod(plot_idx - 1, num_plots_x) == 0
        ylabel('j×B vs grad P [kN/m^3]');
    end
end

% 全体タイトル
sgtitle(['Force Balance Analysis (Radial Component) for Date: ', num2str(date)], 'FontSize', 16, 'FontWeight', 'bold');


% --- 画像保存処理 ---
if isSave
    % 保存先ディレクトリ
    base_dir = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\j×B_calc\\figure';
    date_dir = fullfile(base_dir, num2str(date));
    if ~exist(date_dir, 'dir')
        mkdir(date_dir);
    end

    % ファイル名作成
    if date == 240610
        case_name = 'CaseI';
    elseif date == 240611
        case_name = 'CaseO';
    else
        error('Unsupported date. Please set date to 240610 or 240611.');
    end

    start_str = num2str(start_time_us);
    end_str = num2str(start_time_us + (total_plots-1)*time_step_us);
    save_name = sprintf('ForceBalance_%s_%s_%sus-%sus.png', case_name, num2str(date), start_str, end_str);
    save_path = fullfile(date_dir, save_name);

    saveas(gcf, save_path);
    disp(['Figure saved to: ', save_path]);
end