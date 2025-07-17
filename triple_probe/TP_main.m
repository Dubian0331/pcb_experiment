%% ========================================================================
%  トリプルプローブデータ解析 メイン実行ファイル
% =========================================================================
clear; close all; clc;

%% ------------------------------------------------------------------------
%  ユーザー設定
% -------------------------------------------------------------------------
% --- ケース設定 ---
date = 240610; % case-I
% date = 240611; % case-O

filepath = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\MachProbe_data\';
DOCID = '1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';

if date == 240610
    shotlist = [18, 21, 24, 27, 33, 35, 38, 41]; % case-I (47を除外)
    case_name = 'Case-I';
elseif date == 240611
    shotlist = [46, 49, 52, 55, 58, 61, 66, 69]; % case-O (72を除外)
    case_name = 'Case-O';
else
    error('Unknown date. Please set date to 240610 or 240611.');
end

% --- 解析パラメータ ---
params.V2 = 20;
params.V3 = 10;
params.A = 1.0; % Hガスの原子量
params.S_probe = 1.4e-5; % 電極表面積 [m^2]
params.time_us_unified = 470:0.1:550;
params.time_slice_indices = 4700:5500;
params.lowpass_freq = 1e5; % ローパスフィルタ周波数 [Hz]
params.sampling_freq = 1e7; % データサンプリング周波数 [Hz]
params.smoothing_window = 10; % 電流平滑化のウィンドウサイズ
params.final_smoothing_window = 7; % 最終平滑化のウィンドウサイズ

%% ------------------------------------------------------------------------
%  データ処理の実行
% -------------------------------------------------------------------------
fprintf('Processing data for %s (Date: %d)...\n', case_name, date);

% process_triple_probe_data関数を呼び出し
[Te_matrix_eV, ne_matrix, R_values, time_values, matFileName] = ...
    process_triple_probe_data(date, filepath, DOCID, shotlist, params, case_name);

fprintf('Data processing complete. MAT file saved to: %s\n', matFileName);

%% --- ★★★ 追加: プラズマ圧力の計算とプロット ★★★
fprintf('Calculating and plotting plasma pressure...\n');

% --- 1. 圧力計算の準備 ---
k_B = 1.380649e-23; % ボルツマン定数 [J/K]
K2ev = 11604.5250061657; % eVからKelvinへの変換係数

% 電子温度をeVから絶対温度(Kelvin)に変換
Te_matrix_K = Te_matrix_eV * K2ev;

% --- 2. プラズマ圧力 (電子) を計算 ---
% Pe = ne * kB * Te [単位: Pa]
Pe_matrix = ne_matrix .* (k_B * Te_matrix_K);

% --- 3. 圧力のカラーコンタープロットを作成 ---
figure('Name', [case_name, ' Plasma Pressure Plot']);

% グリッドデータを作成
[T_grid, R_grid] = meshgrid(time_values, R_values);

% データをプロット
contourf(T_grid, R_grid, Pe_matrix, 100, 'LineColor', 'none');

% カラーバーとラベル
cb = colorbar;
cb.Title.String = 'P_e [Pa]';
colormap('jet');
% clim([0, 300]); % 必要に応じてカラーマップの範囲を手動で調整

% 軸ラベルとタイトル
xlabel('Time [\mus]');
ylabel('R [m]');
title([case_name, ' plot of Electron Pressure [Pa]']);
set(gca, 'FontSize', 12);

% 表示範囲 (既存のプロット設定を流用)
xlim([470, 510]);
ylim([0.1, 0.3]);

fprintf('Plasma pressure plot generated.\n');
% --- 追加部分ここまで ---

%% --- ★★★ 追加: 特定時刻での圧力勾配(grad P)のプロット（複数） ★★★

% --- 1. プロットしたい時刻をリストで指定 [µs] ---
target_times = [480, 485, 490, 495]; 

fprintf('指定された各時刻の圧力勾配をプロットします...\n');

% --- 2. 指定した各時刻でループ処理 ---
for i = 1:length(target_times)
    
    current_target_time = target_times(i);
    
    % a) 指定した時刻に最も近い、実際の時間軸上のインデックスを探す
    [~, time_index] = min(abs(time_values - current_target_time));
    actual_time = time_values(time_index); % 実際に抜き出す時刻
    
    % b) その時刻の圧力プロファイル（R方向の分布）を抜き出す
    P_vs_R_profile = Pe_matrix(:, time_index);
    
    % c) R方向の圧力勾配 (dP/dR) を計算
    grad_P = gradient(P_vs_R_profile, R_values);
    
    % d) 新しいfigureを作成してプロット
    figure('Name', ['Pressure Gradient at t = ', num2str(actual_time), ' us']);
    
    plot(R_values, grad_P, '-o', 'LineWidth', 1.5);
    grid on;
    hold on;
    yline(0, 'k--', 'HandleVisibility', 'off'); % ゼロのライン
    hold off;

    % ラベルとタイトル
    xlabel('Radius (R) [m]');
    ylabel('Pressure Gradient (grad P) [Pa/m]');
    title(['Radial Pressure Gradient at t = ', num2str(actual_time), ' \mus']);
    
end

disp('圧力勾配のプロットが完了しました。');
%% ------------------------------------------------------------------------
%  2次元プロットの生成と保存
% -------------------------------------------------------------------------
if any(Te_matrix_eV(:) > 0) % 有効なデータがあるか確認
    fprintf('Generating 2D plots...\n');
    
    % プロット用オプションを設定
    plot_options = struct();
    plot_options.case_name = case_name;
    plot_options.date = date;
    plot_options.Te_clim = [0, 15];
    plot_options.ne_clim = [0, 1.5e20];
    plot_options.time_lim = [470, 510];
    plot_options.R_lim = [0.1, 0.3];
    % plot_options.save_figures = false;
    plot_options.save_path = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\figure';
    
    % グラフを保存するかユーザーに確認するダイアログを表示
    save_question = '生成されたグラフを保存しますか？';
    save_title = '保存確認';
    choice = questdlg(save_question, save_title, 'はい', 'いいえ', 'いいえ');

    % ユーザーの選択に応じて保存フラグを設定
    if strcmp(choice, 'はい')
        plot_options.save_figures = true;
    else
        plot_options.save_figures = false;
    end

    % プロット関数を呼び出し
    plot_triple_probe_data(Te_matrix_eV, ne_matrix, R_values, time_values, plot_options);
    
    fprintf('Plotting complete.\n');
else
    warning('有効なデータが計算されませんでした。プロットをスキップします。');
end