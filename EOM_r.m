%% ========================================================================
%  運動方程式の項の比較 (慣性項 vs JxB力)
% =========================================================================

%% ★★★ ユーザー設定 ★★★
% --- プロットするケースを選択 ---
% Case-I
date = 240610;
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240610023.mat');
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240610\triple_data2D_case-I.mat');
z_index = 20; % プロットするzの位置 (z=0に対応)

% % case-O
% date = 240611;
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240611055.mat');
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240611\triple_data2D_case-O.mat');
% z_index = 22; % プロットするzの位置 (z=0に対応)

% --- プロットのタイミングを設定 ---
start_time_us = 479;
time_step_us = 1;
num_plots_y = 5;
num_plots_x = 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% --- 1. 慣性項 (dvr/dt) の計算 ---
fprintf('Calculating acceleration from magnetic axis motion...\n');

% a) 磁気軸の位置(r)の時系列データを取得
[magAxisList, ~] = get_axis_x_multi(grid2D, data2D);
[~, primary_axis_idx] = max(magAxisList.psi, [], 1);
r_axis_raw = zeros(1, numel(data2D.trange));
for t_idx = 1:numel(data2D.trange)
    if ~isnan(primary_axis_idx(t_idx))
        r_axis_raw(t_idx) = magAxisList.r(primary_axis_idx(t_idx), t_idx);
    else
        r_axis_raw(t_idx) = NaN;
    end
end
time_sec = data2D.trange * 1e-6;

% b) 位置データを平滑化し、速度と加速度を計算
r_axis_smoothed = smoothdata(r_axis_raw, 'movmean', 5, 'omitnan');
vr = gradient(r_axis_smoothed, time_sec); % [m/s]
ar = gradient(vr, time_sec);              % [m/s^2]

% c) 電子密度データを準備
ne_data = squeeze(triple_data2D.ne_avg); % (R, Time)
mi = 1.672e-27; % イオン質量 [kg]

%% --- 2. プロット処理 ---
figure('WindowState', 'maximized');
r_axis_plot = grid2D.rq(:,1); % プロット用のR座標軸
jxb_time_axis = data2D.trange;

total_plots = num_plots_y * num_plots_x;
for plot_idx = 1:total_plots
    time_in_us = start_time_us + (plot_idx - 1) * time_step_us;
    
    subplot(num_plots_y, num_plots_x, plot_idx);
    hold on;
    grid on;
    
    % --- JxB力の計算 ---
    [~, target_time_idx] = min(abs(jxb_time_axis - time_in_us));
    Jr = data2D.Jr(:, :, target_time_idx); Jt = data2D.Jt(:, :, target_time_idx); Jz = data2D.Jz(:, :, target_time_idx);
    Br = data2D.Br(:, :, target_time_idx); Bt = data2D.Bt(:, :, target_time_idx); Bz = data2D.Bz(:, :, target_time_idx);
    J_vectors = [Jr(:), Jt(:), Jz(:)];
    B_vectors = [Br(:), Bt(:), Bz(:)];
    j_B_force_vectors = cross(J_vectors, B_vectors, 2);
    j_B_r_2D = reshape(j_B_force_vectors(:, 1), size(Jr));
    jxb_profile = j_B_r_2D(:, z_index) / 1000; % [kN/m^3]

    jxb_to_plot = (Jt(:, z_index) .* Bz(:, z_index)) / 1000;
    
    % plot(r_axis_plot, jxb_to_plot, 'b-', 'LineWidth', 1.5, 'DisplayName', 'JtBz');
    
    % --- 慣性項の計算 ---
    % a) 現在時刻の加速度を取得
    current_ar = ar(target_time_idx);
    
    % b) 現在時刻の電子密度プロファイルを取得 (補間)
    [~, ne_time_idx] = min(abs(time_values - time_in_us));
    ne_profile_orig = ne_data(:, ne_time_idx);
    ne_profile_aligned = interp1(R_values, ne_profile_orig, r_axis_plot, 'pchip', 'extrap');
    
    % c) 慣性項の力 F_inertia = ni*mi*ar を計算
    inertial_force_profile = (ne_profile_aligned * mi * current_ar) / 1000; % [kN/m^3]
    
    plot(r_axis_plot, inertial_force_profile, 'g--', 'LineWidth', 1.5, 'DisplayName', 'm_in_ia_r');
    
    % --- グラフ装飾 ---
    yline(0, 'k-', 'HandleVisibility', 'off');
    title(['t = ', num2str(time_in_us), ' us']);
    if plot_idx == 1, legend('FontSize', 7); end
    if plot_idx > total_plots - num_plots_x, xlabel('r [m]'); end
    if mod(plot_idx - 1, num_plots_x) == 0, ylabel('Force [kN/m^3]'); end
end
sgtitle(['Inertial Term vs JxB Force for Date: ', num2str(date)], 'FontSize', 16, 'FontWeight', 'bold');