%% ========================================================================
%  圧力勾配(grad P)を Te/ne のグリッドで計算し、MATファイルに保存する
% =========================================================================
clear; close all; clc;
% --- ユーザー設定 ---
% date = 240610;
date = 240611; 
fprintf('Processing date: %d\n', date);

%% --- 1. 元となるTe, neデータをロード ---
if date == 240610 % Case-I
    te_mat_path = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240610\triple_data2D_case-I.mat';
elseif date == 240611 % Case-O
    te_mat_path = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240611\triple_data2D_case-O.mat';
else
    error('Date not recognized.');
end

try
    % ★★★変更点: ここでロードした time_values と R_values を正とする★★★
    load(te_mat_path, 'triple_data2D', 'R_values', 'time_values');
    Te_data_eV = squeeze(triple_data2D.Te_avg);
    ne_data = squeeze(triple_data2D.ne_avg);
    fprintf('Te/ne data loaded successfully.\n');
catch ME
    error('Te/ne .matファイルのロードに失敗しました: %s', ME.message);
end

%% --- 2. 圧力勾配の計算 (補間は不要) ---
% 物理定数
k_B = 1.380649e-23;
K2ev = 11604.5250061657;

% 圧力 Pe = ne * kB * Te を計算
Te_data_K = Te_data_eV * K2ev;
Pe_matrix = ne_data .* (k_B * Te_data_K); % Peの次元は (R, Time)

% 全ての時間について、R方向の圧力勾配を計算
num_times = length(time_values);
grad_P_2D = zeros(size(Pe_matrix)); % (R, Time)

% ★★★変更点: 転置して(Time, R)で計算する方が直感的★★★
Pe_matrix_T = Pe_matrix.'; % (Time, R) に転置
grad_P_2D_T = zeros(size(Pe_matrix_T));

for t_idx = 1:num_times
    P_vs_R_profile = Pe_matrix_T(t_idx, :);
    % R座標ベクトル(R_values)を使って勾配を計算
    grad_P_vs_R = gradient(P_vs_R_profile, R_values);
    grad_P_2D_T(t_idx, :) = grad_P_vs_R;
end

%% --- 3. 計算結果をMATファイルに保存 ---
% 保存する変数を定義 (Te/neの軸をそのまま保存)
time_axis = time_values;
R_axis = R_values;
grad_P_2D = grad_P_2D_T; % (Time, R)の形式で保存

save_folder = fullfile('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat', num2str(date));
if ~exist(save_folder, 'dir'), mkdir(save_folder); end
savename = fullfile(save_folder, ['grad_P_data_', num2str(date), '.mat']);

save(savename, 'grad_P_2D', 'time_axis', 'R_axis');
fprintf('圧力勾配データを保存しました: %s\n', savename);