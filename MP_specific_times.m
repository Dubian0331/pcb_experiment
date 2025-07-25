%% ========================================================================
%  マッハプローブ流速と磁気面のプロット（特定時間指定版）
% =========================================================================
clear;
close all;
clc;

%% ★★★ ユーザー設定 ★★★
% --- 解析対象のケース ---
% Case-I
MP.date = 240610;
disp_shotnum = 18;
MP.shotlist = [18, 20, 23, 26:27, 31:33, 35:38, 40:41, 44, 47];
MP.ng_ch = [3, 4, 5, 7, 8, 9];

% % Case-O
% MP.date = 240611;
% disp_shotnum = 47;
% MP.shotlist = [47, 48, 51:52, 54:60, 62, 65:70];
% MP.ng_ch = [4, 5, 7, 8, 9];

% --- その他の設定 ---
machpathname = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\mach_probe_data';
savepathname = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\mach_probe';

%% --- 1. データ準備 ---
% --- 磁気面データのロード ---
loadfile = fullfile('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data', strcat(num2str(MP.date), sprintf('%03d', disp_shotnum), '.mat'));
load(loadfile, 'data2D', 'grid2D');
disp(strcat('Loaded magnetic data: ', loadfile));

% --- マッハプローブデータ計算のためのパラメータ設定 ---
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';
T=getTS6log(DOCID); T=searchlog(T,'date',MP.date);
MP.rlist=T.MachProbeRPosition_cm_(MP.shotlist)*10;
MP.shot = 3; MP.mesh = 100; MP.meshr = 100; MP.trange = 470:0.1:530;
k = 1.38e-23; q = 1.60e-19; Te = 10 * q / k; mi = 1.672e-27;
MPset.K = 2.5; MPset.cs = sqrt(2 * k * Te / mi);

% --- マッハプローブデータ（流速）の計算 ---
MPdata2D = MP_calc(machpathname, savepathname, MP, MPset);

%% --- 2. プロット関数の呼び出し ---
% プロット設定
plotoptions.times_to_plot = [484, 495, 506];
plotoptions.axisfontsize = 16;
plotoptions.titlefontsize = 18;
plotoptions.isSave = true;

MP_plot_with_psi_specific_times(MP, MPdata2D, data2D, grid2D, savepathname, plotoptions);

%% ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
%  ここより下はローカル関数です。
% ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★

function MP_plot_with_psi_specific_times(MP, MPdata2D, data2D, grid2D, savepath, plotoptions)
    % --- 1. プロット準備 ---
    times_to_plot = plotoptions.times_to_plot;
    num_plots = length(times_to_plot);
    
    % プロット数に応じてレイアウトを自動計算
    max_cols = 4; % 横に並べるプロット数の上限
    if num_plots <= max_cols
        num_rows = 1;
        num_cols = num_plots;
    else
        num_rows = ceil(num_plots / max_cols);
        num_cols = max_cols;
    end

    % Figureサイズをsubplot数に応じて自動調整
    base_width = 400; base_height = 450; top_margin = 100;
    fig_width = base_width * num_cols;
    fig_height = base_height * num_rows + top_margin;
    screen_size = get(0, 'ScreenSize');
    left = (screen_size(3) - fig_width) / 2;
    bottom = (screen_size(4) - fig_height) / 2;

    figure('Position', [left, bottom, fig_width, fig_height], 'visible', 'on');
    tiledlayout(num_rows, num_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

    % --- 2. 指定された時間でループしてプロット ---
    for m = 1:num_plots
        time_in_us = times_to_plot(m);
        
        % 描画したい時刻に最も近いデータインデックスを探す
        [~, idx_mp] = min(abs(MP.trange - time_in_us));
        [~, idx_mag] = min(abs(data2D.trange - time_in_us));
        
        nexttile;
        hold on;
        
        % 流速をカラーコンターでプロット
        contourf(MPdata2D.zq, MPdata2D.rq, squeeze(MPdata2D.flow_forplot(idx_mp, :, :)), 100, 'edgecolor', 'none');
        
        % 磁気面を線コンターでプロット
        contour(grid2D.zq, grid2D.rq, data2D.psi(:, :, idx_mag), -6e-3:0.4e-3:6e-3, 'k', 'LineWidth', 1);
        
        hold off;
        
        % グラフ装飾
        colormap(redblue(300));
        clim([-40 40]);
        title([num2str(time_in_us),' us'], 'FontSize', plotoptions.titlefontsize);
        
        % 軸ラベル (外側のみ表示)
        current_tile = gca;
        if current_tile.Layout.Tile > num_plots - num_cols
            xlabel('z [m]', 'FontSize', plotoptions.axisfontsize);
        end
        if mod(current_tile.Layout.Tile - 1, num_cols) == 0
            ylabel('r [m]', 'FontSize', plotoptions.axisfontsize);
        end
        
        xlim([-0.08 0.08]);
        ylim([0.10 0.30]);
        daspect([1 1 1]); % 縦横比を1:1に
    end
    
    % 全体タイトルと共通カラーバー
    sgtitle(['Flow Velocity & Psi Surfaces for Date: ', num2str(MP.date)], 'FontSize', 20, 'FontWeight', 'bold');
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.FontSize = plotoptions.axisfontsize;
    cb.Label.String = 'Vt [km/s]';
    
    % (保存処理は必要に応じてここに記述)
end