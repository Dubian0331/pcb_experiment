%% ========================================================================
%  圧力勾配(grad P)を全時間について計算し、MATファイルに保存する (軸直接指定版)
% =========================================================================
clear; close all; clc;

%% ★★★ ユーザー設定 (J×Bの計算軸と一致させてください) ★★★

% --- 解析対象の日付 ---
date = 240610; % or 240611

% --- ターゲットとなる時間軸とR軸を直接定義 ---
% JxBの計算で使った 'grid2D.trange' や 'grid2D.rq' と同じ値を設定してください
target_time_axis = (470:1:550); % 例: 400µsから800µsまで1µs刻み
target_R_axis = (0:0.01:0.3);      % 例: 0mから0.3mまで1cm刻み

fprintf('Processing date: %d\n', date);
fprintf('Target Grid: Time [%.1f, %.1f] us, Radius [%.3f, %.3f] m\n', min(target_time_axis), max(target_time_axis), min(target_R_axis), max(target_R_axis));

%% --- 1. 元となるTe, neデータをロード ---
if date == 240610 % Case-I
    te_mat_path = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240610\triple_data2D_case-I.mat';
elseif date == 240611 % Case-O
    te_mat_path = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240611\triple_data2D_case-O.mat';
else
    error('Date not recognized.');
end

try
    load(te_mat_path, 'triple_data2D', 'R_values', 'time_values');
    % 3D配列から2D行列に変換
    Te_data_eV = squeeze(triple_data2D.Te_avg);
    ne_data = squeeze(triple_data2D.ne_avg);
    fprintf('Te/ne data loaded successfully.\n');
catch ME
    error('Te/ne .matファイルのロードに失敗しました: %s', ME.message);
end

%% --- 2. Te/neデータをターゲットグリッドに補間 ---
% 元データ用のグリッドを作成
[T_orig_grid, R_orig_grid] = meshgrid(time_values, R_values);

% 補間先のターゲットグリッドを作成
[T_target_grid, R_target_grid] = meshgrid(target_time_axis, target_R_axis);

% interp2を使ってTe, neデータをターゲットグリッドに補間
fprintf('Interpolating Te and ne data onto the target grid...\n');
Te_aligned_eV = interp2(T_orig_grid, R_orig_grid, Te_data_eV, T_target_grid, R_target_grid, 'linear');
ne_aligned = interp2(T_orig_grid, R_orig_grid, ne_data, T_target_grid, R_target_grid, 'linear');

% 補間後の次元は (R, Time) なので、(Time, R) に転置
Te_aligned_eV = Te_aligned_eV.';
ne_aligned = ne_aligned.';

%% --- 3. 圧力勾配の計算 ---
% 物理定数
k_B = 1.380649e-23;
K2ev = 11604.5250061657;

% 圧力 Pe = ne * kB * Te を計算
Te_aligned_K = Te_aligned_eV * K2ev;
Te_aligned_K(Te_aligned_K < 0) = 0; % 補間による負の値を補正
Pe_matrix = ne_aligned .* (k_B * Te_aligned_K);

% 全ての時間について、R方向の圧力勾配を計算
num_times = length(target_time_axis);
grad_P_2D = zeros(size(Pe_matrix));

for t_idx = 1:num_times
    P_vs_R_profile = Pe_matrix(t_idx, :);
    % R座標ベクトル(target_R_axis)を使って勾配を計算
    grad_P_vs_R = gradient(P_vs_R_profile, target_R_axis);
    grad_P_2D(t_idx, :) = grad_P_vs_R;
end

%% --- 4. 計算結果をMATファイルに保存 ---
% 保存する変数を定義 (ターゲットの軸も一緒に保存)
time_axis = target_time_axis;
R_axis = target_R_axis;

save_folder = fullfile('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat', num2str(date));
if ~exist(save_folder, 'dir'), mkdir(save_folder); end
savename = fullfile(save_folder, ['grad_P_data_', num2str(date), '.mat']);

save(savename, 'grad_P_2D', 'time_axis', 'R_axis');

fprintf('圧力勾配データを保存しました: %s\n', savename);