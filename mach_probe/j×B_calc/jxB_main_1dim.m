clear all

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 複数時間の同時プロット用
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% % Case-I
% date = '240610';
% load('G:\My Drive\lab\lab_data\pre_processed\a039_4290.mat');

% target_times = [91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103];
% target_times = [93, 107, 121];
% target_times = [75, 80, 84, 87, 95, 98, 103];

% target_times = [71, 72, 73, 74, 75, 77, 79, 80, 81];
% target_times = [81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91];
% target_times = [85, 86, 87, 88, 89, 90];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is for presentation figure. It's a caseI.
% you can see the jxB deepest moving bigger side of r.
% target_times = [88, 90, 92, 94, 96, 98, 100];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% case-O
date = '240611';
load('G:\My Drive\lab\lab_data\pre_processed\a039_4374.mat');
target_times = [89, 94, 99];
% % target_times = [82, 84, 86, 88, 90, 92, 94, 96, 98];
% target_times = [75, 79, 83, 87, 91];
% target_times = [93, 94, 95, 96, 97, 98, 99, 100];
% target_times = [100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110];
% target_times = [71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is for presentation figure. It's a caseO.
% target_times = [78, 84, 86, 91, 93, 94, 96];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
% グリッドデータ (r方向)
r = grid2D.rq; % r方向のグリッド


% 各時間ステップごとに計算とプロット
for t_idx = 1:length(target_times)
    target_time = target_times(t_idx);

    % 必要なデータを取得
    Jr = data2D.Jr(:, :, target_time);
    Jt = data2D.Jt(:, :, target_time);
    Jz = data2D.Jz(:, :, target_time);

    Br = data2D.Br(:, :, target_time);
    Bt = data2D.Bt(:, :, target_time);
    Bz = data2D.Bz(:, :, target_time);

    % 3次元ベクトルを構築
    J = cat(3, Jr, Jt, Jz); % (20×40×3)
    B = cat(3, Br, Bt, Bz); % (20×40×3)

    % 外積計算
    j_B_force = zeros(size(J)); % 外積結果を格納
    for i = 1:size(J, 1)
        for j = 1:size(J, 2)
            j_B_force(i, j, :) = cross(squeeze(J(i, j, :)), squeeze(B(i, j, :))) / 1000;
        end
    end

    % r成分を抽出 (2Dデータ)
    j_B_r = j_B_force(:, :, 1); % r成分

    % 径方向成分
    % プロット (z = 0) caseI
    % plot(r(:,1), j_B_r(:, 20), 'DisplayName', ['t=', num2str(target_time+399), 'us'], 'LineWidth', 3); % 凡例用ラベル
    % % プロット (z = 0) caseO
    % plot(r(:,1), j_B_r(:, 22), 'DisplayName', ['t=', num2str(target_time+399), 'us'], 'LineWidth', 3); % 凡例用ラベル

    j_B_t = j_B_force(:, :, 2); % toroidal成分
    % トロイダル方向成分
    % caseI
    % plot(r(:,1), j_B_t(:, 20), 'DisplayName', ['t=', num2str(target_time+399), 'us'], 'LineWidth', 3); % 凡例用ラベル
    % プロット (z = 0) caseO
    plot(r(:,1), j_B_t(:, 22), 'DisplayName', ['t=', num2str(target_time+399), 'us'], 'LineWidth', 3); % 凡例用ラベル

    legend('FontSize', 10, 'Location', 'southwest')
    hold on

end

% グラフ
% title('case - I', 'FontSize', 14);
xlabel('r [m]', 'FontSize', 14);
% ylabel('j_{t}B_{z} [kN/m^{3}]', 'FontSize', 14);
yline(0, '--', 'HandleVisibility', 'off'); % y=0ライン
ylim([-35 20])
% % legend('time=488us', 'time=493us', 'time=498us', 'Location', 'best'); % 凡例を自動配置


% 図を保存
saveFolder = "G:\My Drive\lab\lab_data\mach_probe\j×B_calc\figure";
saveas(gcf, strcat(saveFolder, '\', num2str(date), '\jt_Bz_force-r_multiple_times_caseO_toroidal', num2str(min(target_time)+399), 'us-', num2str(max(target_time)+399), 'us'), 'png')

fprintf("save your file %s\n", num2str(date));

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 同じウィンドウにsubplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% Case-I
% date = '240610';
% load('G:\My Drive\lab\lab_data\pre_processed\a039_4290.mat');
% target_times = [79, 80, 81, 82, 83, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103];
% target_times = [86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110];

% case-O
date = '240611';
load('G:\My Drive\lab\lab_data\pre_processed\a039_4374.mat');
% target_times = [89, 94, 99];
target_times = [87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100];


% グリッドデータ (r方向)
r = grid2D.rq; % r方向
z = grid2D.zq;

% subplot 行数と列数の計算
num_plots = length(target_times);
num_rows = ceil(sqrt(num_plots));
num_cols = ceil(num_plots / num_rows);

% プロット用のウィンドウを作成
figure;
set(gcf, 'Position', get(0, 'Screensize')); % ウィンドウを最大化

% 各時間ステップごとに計算とプロット
for t_idx = 1:num_plots
    target_time = target_times(t_idx);

    % 必要なデータを取得
    Jr = data2D.Jr(:, :, target_time);
    Jt = data2D.Jt(:, :, target_time);
    Jz = data2D.Jz(:, :, target_time);

    Br = data2D.Br(:, :, target_time);
    Bt = data2D.Bt(:, :, target_time);
    Bz = data2D.Bz(:, :, target_time);

    % 3次元ベクトルを構築
    J = cat(3, Jr, Jt, Jz); % (20×40×3)
    B = cat(3, Br, Bt, Bz); % (20×40×3)

    % 外積計算
    j_B_force = zeros(size(J)); % 外積結果を格納
    for i = 1:size(J, 1)
        for j = 1:size(J, 2)
            j_B_force(i, j, :) = cross(squeeze(J(i, j, :)), squeeze(B(i, j, :))) / 1000;
        end
    end

    % r成分を抽出 (2Dデータ)
    jt_Bz_r = j_B_force(:, :, 1); % r成分

    % subplot を作成してプロット
    subplot(num_rows, num_cols, t_idx);
    plot(r(:,1), jt_Bz_r(:, 22), 'DisplayName', ['t=', num2str(target_time+399), 'us']);
    title(['t = ', num2str(target_time+399), ' us']);
    xlabel('r [m]');
    xlim([0.1 0.3])
    ylabel('j_{t}B_{z} [kN/m^{3}]');
    ylim([-40 30])
    yline(0, 'HandleVisibility', 'off'); % y=0ライン
    % legend('Location', 'northwest')
end

% 全体のタイトルを追加
sgtitle('case - O');

% % 図を保存
saveFolder = "G:\My Drive\lab\lab_data\mach_probe\j×B_calc\figure";
saveas(gcf, strcat(saveFolder, '\', num2str(date), '\jt_Bz_force-r_subplots_caseO_', num2str(min(target_times)+399), 'us-', num2str(max(target_times)+399), 'us'), 'png');

fprintf("save your file %s\n", num2str(date));

%%
% 3成分を個別にプロット
figure;

% r方向成分
subplot(1, 3, 1);
surf(r, z, j_B_force(:, :, 1)); % r成分
shading interp;
colorbar;
title('jxB_{r}');
xlabel('r [m]');
ylabel('z [m]');
zlabel('jxB_{r} [kN/m^{3}]');

% z方向成分
subplot(1, 3, 2);
surf(r, z, j_B_force(:, :, 2)); % z成分
shading interp;
colorbar;
title('jxB_{z}');
xlabel('r [m]');
ylabel('z [m]');
zlabel('jxB_{z} [kN/m^{3}]');

% t方向成分
subplot(1, 3, 3);
surf(r, z, j_B_force(:, :, 3)); % t成分
shading interp;
colorbar;
title('j_{B,t}');
xlabel('r [m]');
ylabel('z [m]');
zlabel('j_{B,t} [kN/m^{3}]');

% ベクトルの大きさを計算
j_B_magnitude = sqrt(j_B_force(:, :, 1).^2 + j_B_force(:, :, 2).^2 + j_B_force(:, :, 3).^2);

% 3次元プロット
figure;
surf(r, z, j_B_magnitude); % ベクトルの大きさをプロット
shading interp; % 補間で滑らかに表示
colorbar; % カラーバーを追加
title('|jxB| Magnitude');
xlabel('r [m]');
ylabel('z [m]');
zlabel('|jxB| [kN/m^{3}]');

% 保存
% saveas(gcf, strcat(saveFolder, '\', num2str(date), '\j_B_force_components_t', num2str(target_time+399), 'us'), 'png');
% fprintf("3D component plots saved for date %s\n", num2str(date));
