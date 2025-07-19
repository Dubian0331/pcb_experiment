%% ========================================================================
%  トリプルプローブデータ解析 メイン実行ファイル
% =========================================================================
% clear; close all; clc;

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