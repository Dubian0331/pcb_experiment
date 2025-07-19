% clear all;
% 
% % データのロード
% date = '240610';
% load('G:\My Drive\lab\lab_data\pre_processed\a039_4290.mat');
% 
% % 対象時間
% target_times = [88, 90, 92, 94, 96, 98, 100];
% 
% % グリッドデータ (r, z)
% r = grid2D.rq; % r方向のグリッド
% z = grid2D.zq; % z方向のグリッド
% 
% % フィギュア作成
% figure;
% num_plots = length(target_times); % サブプロット数
% 
% for t_idx = 1:num_plots
%     % 現在の時間を取得
%     target_time = target_times(t_idx);
% 
%     % 必要なデータを取得
%     Jr = data2D.Jr(:, :, target_time);
%     Jt = data2D.Jt(:, :, target_time);
%     Jz = data2D.Jz(:, :, target_time);
% 
%     Br = data2D.Br(:, :, target_time);
%     Bt = data2D.Bt(:, :, target_time);
%     Bz = data2D.Bz(:, :, target_time);
% 
%     % 3次元ベクトルを構築
%     J = cat(3, Jr, Jt, Jz); % (20×40×3)
%     B = cat(3, Br, Bt, Bz); % (20×40×3)
% 
%     % 外積計算
%     j_B_force = zeros(size(J)); % 外積結果を格納
%     for i = 1:size(J, 1)
%         for j = 1:size(J, 2)
%             j_B_force(i, j, :) = cross(squeeze(J(i, j, :)), squeeze(B(i, j, :))) / 1000;
%         end
%     end
% 
%     % r成分を抽出 (2Dデータ)
%     jt_Bz_r = j_B_force(:, :, 1); % r成分
% 
%     % 等高線レベルを自動計算
%     % contour_levels = linspace(min(jt_Bz_r(:)), max(jt_Bz_r(:)), 10);
% 
%     % 等高線レベルをユーザーに出力
%     % fprintf('Time %d us - Contour Levels: [%s]\n', target_time + 399, num2str(contour_levels));
% 
%     % サブプロット
%     subplot(ceil(sqrt(num_plots)), ceil(sqrt(num_plots)), t_idx);
%     hold on
%     % contourf(z, r, jt_Bz_r, contour_levels, 'LineStyle', 'none');
%     contourf(z, r, jt_Bz_r, 100, 'LineStyle','none');
%     contourf(grid2D.zq, grid2D.rq, data2D.psi(:, :, target_time), -6e-3:0.4e-3:6e-3, 'k', 'Fill', 'off', 'LineWidth', 1);
%     colorbar;
%     title(['t = ', num2str(target_time + 399), ' us']);
%     xlabel('z [m]');
%     ylabel('r [m]');
%     colormap("whitejet")
%     clim([-60, 60]);
%     % clim([min(jt_Bz_r(:)), max(jt_Bz_r(:))]); % カラーバーの範囲
% end



%%
clear all;

% データのロード
% case-I
% date = '240610';
% load('G:\My Drive\lab\lab_data\pre_processed\a039_4290.mat');
% target_times = [93];
% target_times = [88, 90, 92, 94, 96, 98];
% target_times = [93, 94, 107, 121];


% % case-O
date = '240611';
load('G:\My Drive\lab\lab_data\pre_processed\a039_4374.mat');
target_times = [89];
% target_times = [78, 84, 86, 91, 93, 94, 96];
% target_times = [78, 84, 86, 91, 93, 94];
% target_times = [78, 91, 94];
% target_times = [88, 89, 90, 99];


% グリッドデータ (r, z)
r = grid2D.rq; % r方向のグリッド
z = grid2D.zq; % z方向のグリッド

% フィギュア作成
figure;
set(gcf, 'Position', get(0, 'Screensize')); % ウィンドウを最大化
num_plots = length(target_times); % サブプロット数
n_cols = ceil(sqrt(num_plots)); % 横の列数
n_rows = ceil(num_plots / n_cols); % 縦の行数

% データ範囲の初期化（全体カラーバーの範囲を統一）
global_min = Inf;
global_max = -Inf;

% 等高線レベルの計算のために全データの最小・最大を算出
for t_idx = 1:num_plots
    target_time = target_times(t_idx);
    Jr = data2D.Jr(:, :, target_time);
    Jt = data2D.Jt(:, :, target_time);
    Jz = data2D.Jz(:, :, target_time);
    Br = data2D.Br(:, :, target_time);
    Bt = data2D.Bt(:, :, target_time);
    Bz = data2D.Bz(:, :, target_time);
    J = cat(3, Jr, Jt, Jz);
    B = cat(3, Br, Bt, Bz);
    j_B_force = zeros(size(J));
    for i = 1:size(J, 1)
        for j = 1:size(J, 2)
            j_B_force(i, j, :) = cross(squeeze(J(i, j, :)), squeeze(B(i, j, :))) / 1000;
        end
    end

    jt_Bz_r = j_B_force(:, :, 1);
    j_B_toroidal = j_B_force(:, :, 2);
    global_min = min(global_min, min(jt_Bz_r(:)));
    global_max = max(global_max, max(jt_Bz_r(:)));
end

% 等高線レベルを計算
contour_levels = linspace(global_min, global_max, 10);

% サブプロットのプロット処理
for t_idx = 1:num_plots
    target_time = target_times(t_idx);

    Jr = data2D.Jr(:, :, target_time);
    Jt = data2D.Jt(:, :, target_time);
    Jz = data2D.Jz(:, :, target_time);
    Br = data2D.Br(:, :, target_time);
    Bt = data2D.Bt(:, :, target_time);
    Bz = data2D.Bz(:, :, target_time);
    J = cat(3, Jr, Jt, Jz);
    B = cat(3, Br, Bt, Bz);
    j_B_force = zeros(size(J));
    for i = 1:size(J, 1)
        for j = 1:size(J, 2)
            j_B_force(i, j, :) = cross(squeeze(J(i, j, :)), squeeze(B(i, j, :))) / 1000;
        end
    end
    jt_Bz_r = j_B_force(:, :, 1);
    j_B_toroidal = j_B_force(:, :, 2);

    % サブプロットの作成
    subplot(n_rows, n_cols, t_idx);
    hold on
    % contourf(z, r, jt_Bz_r, 100, 'LineStyle', 'none');
    % contourf(z(17:34), r, j_B_toroidal(:, 17:34), 100, 'LineStyle', 'none');
    contourf(z(:, 17:24), r(:, 17:24), j_B_toroidal(:, 17:24), 100, 'LineStyle', 'none');
    contourf(grid2D.zq, grid2D.rq, data2D.psi(:, :, target_time), -20e-3:0.4e-3:60e-3, 'k', 'Fill', 'off', 'LineWidth', 1);    

    % 軸メモリとラベルの設定
    if mod(t_idx - 1, n_cols) == 0
        ylabel('r [m]', 'FontSize', 20); % 一番左の列のみ縦軸を表示
    else
        set(gca, 'YTickLabel', []); % その他は非表示
    end

    if t_idx > (n_rows - 1) * n_cols
        xlabel('z [m]', 'FontSize', 20); % 一番下の行のみ横軸を表示
    else
        set(gca, 'XTickLabel', []); % その他は非表示
    end

    pbaspect([2 1 1]);
    set(gca, 'FontSize', 15);
    title(['t = ', num2str(target_time + 399), ' us'], 'FontSize', 20);
end

% 全体のカラーバーを作成
h = colorbar('Position', [0.93 0.1 0.02 0.8]); % フィギュアの右側に配置
set(h, 'FontSize', 15);
ylabel(h, 'j×B_{toroidal} [kN/m^{3}]', 'FontSize', 20);
clim([-40, 40]); % 全体カラーバーの範囲を設定
colormap("whitejet")


% 図を保存
saveFolder = "G:\My Drive\lab\lab_data\mach_probe\j×B_calc\figure";
saveas(gcf, strcat(saveFolder, '\', num2str(date), '\dim2_toroidal', num2str(min(target_times)+399), 'us-', num2str(max(target_times)+399), 'us'), 'png')

fprintf("save your file %s\n", num2str(date));