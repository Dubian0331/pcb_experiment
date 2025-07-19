%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 磁気軸のR方向の移動からプラズマのR方向への移動速度を見積もる
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
% --- プロットするケースを選択 ---
% Case-I
date = 240610;
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240610023.mat');
z_index = 20; % プロットするzの位置 (z=0に対応)

% % case-O
% date = 240611;
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240611055.mat');
% z_index = 22; % プロットするzの位置 (z=0に対応)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time_us = 480; % プロットを開始する時刻 [us]
time_step_us = 1;    % プロットの時間間隔 [us]
num_plots_y = 5;     % 縦に並べるプロット数
num_plots_x = 6;     % 横に並べるプロット数
dispGradP = false; % gradPを表示するか
dispSum = false; % J×BとgradPの合計値を表示するか
isSave = true; % 画像を保存するかどうか
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 
% --- ステップ1: 磁気軸の位置(r)の時系列データを取得 ---
[magAxisList, xPointList] = get_axis_x_multi(grid2D, data2D);

% 2つある磁気軸のうち、psiが大きい方（主磁気軸）を選択
[~, primary_axis_idx] = max(magAxisList.psi, [], 1);
% 各時刻(列)について、主磁気軸のインデックスを使ってr座標を取得
r_axis_raw = zeros(1, numel(data2D.trange));
for t_idx = 1:numel(data2D.trange)
    if ~isnan(primary_axis_idx(t_idx))
        r_axis_raw(t_idx) = magAxisList.r(primary_axis_idx(t_idx), t_idx);
    else
        r_axis_raw(t_idx) = NaN;
    end
end
% 時間軸をマイクロ秒から秒に変換
time_sec = data2D.trange * 1e-6;

% --- ステップ2: 位置データを平滑化 ---
% ノイズの影響を抑えるため、微分前にスムージングを行う（ウィンドウサイズは要調整）
r_axis_smoothed = smoothdata(r_axis_raw, 'movmean', 5, 'omitnan');

% --- ステップ3: 移動速度 (vr) を計算 [m/s] ---
vr = gradient(r_axis_smoothed, time_sec);

% --- ステップ4: 加速度 (ar) を計算 [m/s^2] ---
ar = gradient(vr, time_sec);

%% --- 結果のプロット ---
figure('WindowState', 'maximized');

% プロット範囲（470us～530us）に制限
plot_start_us = 500;
plot_end_us = 520;
plot_mask = (time_sec * 1e6 >= plot_start_us) & (time_sec * 1e6 <= plot_end_us);

% 1. 位置のプロット
subplot(3, 1, 1);
plot(time_sec(plot_mask) * 1e6, r_axis_raw(plot_mask), 'color', [0.7 0.7 0.7], 'DisplayName', 'Raw Position');
hold on;
plot(time_sec(plot_mask) * 1e6, r_axis_smoothed(plot_mask), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Smoothed Position');
grid on;
ylabel('Magnetic Axis r-pos [m]');
title('Magnetic Axis Position');
legend;

% 2. 速度のプロット（y軸-20～20に設定）
subplot(3, 1, 2);
plot(time_sec(plot_mask) * 1e6, vr(plot_mask) / 1000, 'g-', 'LineWidth', 1.5); % 単位を km/s に変換
grid on;
ylabel('Radial Velocity (v_r) [km/s]');
title('Estimated Radial Velocity');

% 3. 加速度のプロット
subplot(3, 1, 3);
plot(time_sec(plot_mask) * 1e6, ar(plot_mask), 'r-', 'LineWidth', 1.5);
grid on;
ylabel('Radial Acceleration (a_r) [m/s^2]');
title('Estimated Radial Acceleration');
xlabel('Time [\mus]');

sgtitle('Estimation of Ion Acceleration from Magnetic Axis Motion', 'FontSize', 14, 'FontWeight', 'bold');