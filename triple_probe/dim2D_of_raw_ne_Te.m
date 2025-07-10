%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case-I, Case-Oのトリプルプローブデータから
% 2次元のTeとneを計算し、プロットする
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

% ファイルのパス設定
% date = 240610; % case-I
date = 240611; % case-O
% filepath = "G:\My Drive\lab\lab_data\mach_probe_rawdata\"; % ファイルパス
filepath = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data\\';
% shotlist = [18, 21, 24, 27, 33, 35, 38, 41, 47]; % case-I
% shotlist = [18, 21, 24, 27, 33, 35, 38, 41]; % case-I
% shotlist = [46, 49, 52, 55, 58, 61, 66, 69, 72]; % case-O
shotlist = [46, 49, 52, 55, 58, 61, 66, 69]; % case-O
% shotlist = [72];
% スプレッドシートからR方向データを取得
DOCID = '1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw'; % スプレッドシートのID
T = getTS6log(DOCID);
node = 'date';
T = searchlog(T, node, date);

% ショット番号に対応するインデックスを取得
shot_indices = arrayfun(@(x) find(T.shot == x, 1), shotlist);

% r_list を取得し、サイズを shotlist と一致させる
r_list = T.tripleProbeRPosition_cm_(shot_indices) * 1e-2; % R座標[m]


% パラメータ設定
V2 = 20; % 電位設定
V3 = 10;
A = 1.0;  % Hガスの原子量
% A = 39.95; % Arガスの原子量
% S_probe = 2.89e-5; % 電極表面積 [m^2]
S_probe = 1.4e-5;
q = 1.60217663e-19; % 電子電荷
kb = 1.380649e-23; % ボルツマン定数
K2ev = 11604.5250061657; % kelvin to eV
mi = A * 1.66054e-27; % イオン質量


% 時間データの統一と行列初期化
num_shots = length(shotlist);
time_us = linspace(470, 550, 300); % 統一された時間軸 (450-600 μs, 300点)
Te_matrix = zeros(num_shots, length(time_us)); % Teデータ格納
ne_matrix = zeros(num_shots, length(time_us)); % neデータ格納

% シンボリック変数宣言
syms I1 I2 I3 Te ne

% Te解法
eqn_Te = (I1 + I2) / (I1 + I3) == ...
    (1 - exp(-q * V2 / (kb * Te))) / (1 - exp(-q * V3 / (kb * Te)));
Te_symbolic = solve(eqn_Te, Te);

% ne解法
eqn_ne = exp(-0.5) * S_probe * ne * q * sqrt(kb * Te / mi) == ...
    (I3 - I2 * exp(-q * (V3 - V2) / (kb * Te))) / (1 - exp(-q * (V3 - V2) / (kb * Te)));
ne_symbolic = solve(eqn_ne, ne);

% ショットごとのデータ読み込みと計算
valid_idx = 0; % 有効なショット数カウンタ
for idx = 1:num_shots
    shot = shotlist(idx);
    filename = strcat(filepath, num2str(date), '\ES_', num2str(date), sprintf('%03d', shot), '.csv');
    
    if exist(filename, 'file')
        valid_idx = valid_idx + 1; % 有効なデータのカウンタを増やす
        % データ読み込み
        data = readmatrix(filename);

        % プロット用
        time_raw = data(:, 1); % ファイル全体の時間データ (プロット用)
        I2_raw_values = data(:, 36); % ファイル全体の列36のデータ (プロット用)
        I3_raw_values = data(:, 37); % ファイル全体の列37のデータ (プロット用)

        time = data(4700:5500, 1); % 時間データ
        I2_values = data(4700:5500, 36); % 列36のデータ
        I3_values = data(4700:5500, 37); % 列37のデータ

        % 不適切なデータの削除
        current_noise_threshold = 1e-6; % 例: 1 microAmpere (A)以下はノイズと見なす
        I2_values(I2_values < 0 | abs(I2_values) < current_noise_threshold) = NaN;
        I3_values(I3_values < 0 | abs(I3_values) < current_noise_threshold) = NaN;
        
        % 欠損データの補完とスムージング (Te/ne計算前の電流データは維持)
        % ここはそのまま残すことで、Te/ne計算の入力電流をクリーンに保つ
        I2_values = fillmissing(I2_values, 'linear');
        I3_values = fillmissing(I3_values, 'linear');
        I2_values = lowpass(I2_values, 0.1e6, 1e7);
        I3_values = lowpass(I3_values, 0.1e6, 1e7);
        I2_values = smoothdata(I2_values, 'movmedian', 10);
        I3_values = smoothdata(I3_values, 'movmedian', 10);
        I1_values = I2_values + I3_values;
        Denominator_Left_Side = I2_values + 2 * I3_values;

       % Te 計算
        try
            % 電流値のコピーを作成し、それに対して操作を行う
            I1_for_calc = I1_values;
            I2_for_calc = I2_values;
            I3_for_calc = I3_values;
        
            % ゼロ除算を避けるための閾値 (非常に小さい値)
            denominator_threshold = 1e-12; 
            
            % --- 新たな問題点: log の引数内の分母 (I2 - I3) のチェック ---
            % I2 - I3 が閾値以下になる（ゼロに近い）インデックスを見つける
            problematic_indices_I2_minus_I3 = find(abs(I2_for_calc - I3_for_calc) < denominator_threshold);
            if ~isempty(problematic_indices_I2_minus_I3)
                I1_for_calc(problematic_indices_I2_minus_I3) = NaN;
                I2_for_calc(problematic_indices_I2_minus_I3) = NaN;
                I3_for_calc(problematic_indices_I2_minus_I3) = NaN;
            end
        
            % --- 新たな問題点: log の引数 (I1 + I3) / (I2 - I3) が非正になる場合のチェック ---
            % まず、I2 - I3 が NaN でない部分の比率を計算
            ratio_arg_for_log = (I1_for_calc + I3_for_calc) ./ (I2_for_calc - I3_for_calc);
            
            % 比率がゼロ以下になるインデックスを見つける (log(0) や log(負の値) を避けるため)
            problematic_indices_log_arg = find(ratio_arg_for_log <= 0);
            if ~isempty(problematic_indices_log_arg)
                I1_for_calc(problematic_indices_log_arg) = NaN;
                I2_for_calc(problematic_indices_log_arg) = NaN;
                I3_for_calc(problematic_indices_log_arg) = NaN;
            end
            
            % --- 新たな問題点: Te_symbolic の最終分母がゼロになる場合 (log の結果が 0) のチェック ---
            % これは (I1 + I3) / (I2 - I3) が 1 になる場合。
            % 非常にまれなケースですが、厳密には以下のようにチェックできます。
            % ただし、浮動小数点数での比較なので、完全に1.0ではなく、非常に1.0に近い範囲をNaNとする
            problematic_indices_final_log_zero = find(abs(ratio_arg_for_log - 1) < denominator_threshold);
            if ~isempty(problematic_indices_final_log_zero)
                I1_for_calc(problematic_indices_final_log_zero) = NaN;
                I2_for_calc(problematic_indices_final_log_zero) = NaN;
                I3_for_calc(problematic_indices_final_log_zero) = NaN;
            end
        
            % Te の計算
            Te_values_raw_from_solve = real(double(subs(Te_symbolic, {I1, I2, I3}, {I1_for_calc, I2_for_calc, I3_for_calc})));

            Te_values = Te_values_raw_from_solve;
            max_physical_Te = 100 * K2ev; % 100 eV 以上は異常値とみなす。
            Te_values(Te_values < 0 | Te_values > max_physical_Te) = NaN; 
        
        catch ME 
            disp(['An error occurred during Te calculation for shot ', num2str(shot), ':']);
            disp(ME.message);
            Te_values = NaN(size(I1_values));
        end

        % ne 計算
        try
            ne_values = real(double(subs(ne_symbolic, {Te, I2, I3}, {Te_values, I2_values, I3_values})));
            ne_values(ne_values < 0) = NaN;
        catch
            % ne_values = zeros(size(I2_values));
            ne_values = NaN(size(I2_values));
        end

        % Te_values(isnan(Te_values)) = 0.1 * K2ev;
        % ne_values(isnan(ne_values)) = 1e17;

        % 補間して時間軸を統一 (NaNを0ではなくNaNで外挿)
        Te_matrix(valid_idx, :) = interp1(time, Te_values, time_us, 'linear', 0); % ここを変更
        ne_matrix(valid_idx, :) = interp1(time, ne_values, time_us, 'linear', 0); % ここを変更

        Te_matrix(valid_idx, :) = smoothdata(Te_matrix(valid_idx, :), 'movmedian', 7); % ウィンドウサイズは調整
        ne_matrix(valid_idx, :) = smoothdata(ne_matrix(valid_idx, :), 'movmedian', 7); % ウィンドウサイズは調整
         
        if shot == shotlist(4) % あるいは確認したい特定のショット番号を指定 (例: if shot == 46)
            figure('Name', ['Raw I2 and I3 for Shot ', num2str(shot)]);
            plot(time_raw, I2_raw_values, 'b-', 'DisplayName', 'I2 (Raw)');
            hold on;
            plot(time_raw, I3_raw_values, 'r-', 'DisplayName', 'I3 (Raw)');
            grid on;
            % 補間後のデータ (Te, ne計算に使用された範囲)
            plot(time, I2_values, 'c-', 'LineWidth', 1.5, 'DisplayName', 'I2 (Corrected)');
            plot(time, I3_values, 'm-', 'LineWidth', 1.5, 'DisplayName', 'I3 (Corrected)');

            xlim([450, 520]);
            xlabel('Time [μs]');
            ylabel('Current [A]');
            title(['Raw I2 and I3 for Shot ', num2str(shot), ' (R = ', num2str(r_list(idx)), ' m)']);
            legend show;
            hold off;
            % プロット後、実行を一時停止して確認したい場合は以下のコメントを外す
            % pause; % Enterキーを押すと次のショットに進む
            CurrentRatio = (I1_values + I2_values) ./ (I1_values + I3_values);
            figure; plot(time, CurrentRatio); title(['Current Ratio for Shot ', num2str(shot)]);
            xlim([450, 520]);
            ylim([-5, 5]); % 適宜範囲を調整
            grid on;
        end
    else
        warning('File %s not found. Skipping this shot.', filename);
    end
end
%%

% 有効データのみを使用
Te_matrix = Te_matrix(1:valid_idx, :);
ne_matrix = ne_matrix(1:valid_idx, :);
r_list = r_list(1:valid_idx);

% 0をNaNに置き換え (Te_valuesが0になることは通常物理的に意味がない場合)
% Te_matrix(Te_matrix == 0) = NaN;
% ne_matrix(ne_matrix == 0) = NaN;


% --------------------------------------------------------------------------------------
% ここから以下の fillmissing はコメントアウトまたは削除し、NaNを保持する
% --------------------------------------------------------------------------------------

% NaN値の補間 (時間方向) -- ここをコメントアウト
for i_shot = 1:size(Te_matrix, 1)
    Te_matrix(i_shot, :) = fillmissing(Te_matrix(i_shot, :), 'linear'); 
    ne_matrix(i_shot, :) = fillmissing(ne_matrix(i_shot, :), 'linear');
end

% NaN値を補完 (線形補完) -- ここもコメントアウトし、NaNをそのまま渡す
% Te_matrix_fill = fillmissing(Te_matrix, 'linear', 1); % R方向で補完
% ne_matrix_fill = fillmissing(ne_matrix, 'linear', 1); % R方向で補完

% 上記の fillmissing をしない場合、Te_matrix_fill は Te_matrix そのものになる
Te_matrix_fill = Te_matrix; % <--- ここを変更
ne_matrix_fill = ne_matrix; % <--- ここを変更


% ----------------------------
% 同じr_listの値でTeとneの平均値を計算
% ----------------------------
unique_r = unique(r_list); % 一意なR値
Te_avg_matrix = zeros(length(unique_r), length(time_us));
ne_avg_matrix = zeros(length(unique_r), length(time_us));

for i = 1:length(unique_r)
    idx_r = (r_list == unique_r(i)); % 現在のR値に対応するインデックス
    % NaN は 'omitnan' で無視されるが、全てNaNの場合は結果もNaNになる
    Te_avg_matrix(i, :) = mean(Te_matrix_fill(idx_r, :), 1, 'omitnan'); % 平均
    ne_avg_matrix(i, :) = mean(ne_matrix_fill(idx_r, :), 1, 'omitnan'); % 平均
end

% --------------------------------------------------------------------------------------
% 平均化後の最終スムージング (オプションだが、試す価値あり) -- ここもコメントアウト
% --------------------------------------------------------------------------------------
% % 時間方向のスムージング (各R位置について)
% for i = 1:size(Te_avg_matrix, 1)
%     Te_avg_matrix(i, :) = smoothdata(Te_avg_matrix(i, :), 'movmedian', 7); % ウィンドウサイズは調整
%     ne_avg_matrix(i, :) = smoothdata(ne_avg_matrix(i, :), 'movmedian', 7); % ウィンドウサイズは調整
% end
% 
% % R方向のスムージング (各時間について)
% % dim = 1 を指定してR方向（行方向）に適用
% Te_avg_matrix = smoothdata(Te_avg_matrix, 'movmedian', 3, 'omitnan', 1);
% ne_avg_matrix = smoothdata(ne_avg_matrix, 'movmedian', 3, 'omitnan', 1); 


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matファイルの作製
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 必要なデータを準備
R_values = unique_r;             % R方向の値 (10個)
time_values = time_us;           % 時間軸 (300点)
Te_data_2D = Te_avg_matrix ./ K2ev; % Teデータを[eV]に変換した2次元データ (10×300)
ne_data_2D = ne_avg_matrix;       % neデータ (10×300)

% 3次元データ構造を作成
triple_data2D = struct(); % 空のストラクトを作成

% Teデータを格納 (サイズ: 10×1×300)
triple_data2D.Te = zeros(length(R_values), 1, length(time_values));
triple_data2D.Te(:, 1, :) = reshape(Te_data_2D, [length(R_values), 1, length(time_values)]);

% neデータを格納 (サイズ: 10×2×300)
triple_data2D.ne = zeros(length(R_values), 1, length(time_values));
triple_data2D.ne(:, 1, :) = reshape(ne_data_2D, [length(R_values), 1, length(time_values)]);


% MATファイルとして保存
targetFolderPath = fullfile('C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\\triple_probe\\mat\\',num2str(date));

if ~exist(targetFolderPath, 'dir')
    mkdir(targetFolderPath);
end

% ファイル名の動的決定
if date == 240610
    matFileName = fullfile(targetFolderPath, 'triple_data2D_case-I.mat');
elseif date == 240611
    matFileName = fullfile(targetFolderPath, 'triple_data2D_case-O.mat');
else
    matFileName = fullfile(targetFolderPath, 'triple_data2D_unknown_case.mat');
end

% 指定したパスにMATファイルを保存
save(matFileName, 'triple_data2D', 'R_values', 'time_values');

% 保存されたファイルのパスを表示
disp(['MATファイル ', matFileName, ' を保存しました。']);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conturfのプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATファイルを読み込む
load(matFileName); % triple_data2D, R_values, time_values を読み込む

% 必要なデータを取り出す
Te_data_2D = squeeze(triple_data2D.Te(:, 1, :)); % 3次元データから2次元データに変換 (10×300)
ne_data_2D = squeeze(triple_data2D.ne(:, 1, :)); % 3次元データから2次元データに変換 (10×300)


% % グリッドデータを作成
[T_grid, R_grid] = meshgrid(time_values, R_values); % 時間 (X軸) と R方向 (Y軸) のグリッド

% % Teプロット作成
figure;
contourf(T_grid, R_grid, Te_data_2D, 10000, 'LineColor', 'none');
clim([0, 10]); % 必要に応じてNaNの範囲も表示されるように調整
xlim([475, 520]);
% カラーバーを追加してタイトルを設定
cb = colorbar;
cb.Title.String = 'T_{e} [eV]';
% 個別に文字サイズを設定
xlabel('Time [\mus]', 'FontSize', 12);     % X軸ラベル
ylabel('R [m]', 'FontSize', 14);           % Y軸ラベル
title('Case - I plot of Te [eV]', 'FontSize', 16); % グラフタイトル
cb.Title.FontSize = 10;                    % カラーバーのタイトル
cb.FontSize = 10;                          % カラーバーの目盛りフォントサイズ
% 軸目盛りフォントサイズを設定
ax = gca; % 現在の座標軸を取得
ax.FontSize = 12; % 軸目盛りフォントサイズを設定
% カラーマップ設定
colormap('jet');
ylim([0.13, 0.3]);
xlim([470, 510]); 
% 保存する
% saveDir = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\figure';
% dateFolder = fullfile(saveDir, num2str(date));
% if ~exist(dateFolder, 'dir')
%     mkdir(dateFolder);
% end
% saveas(gcf, fullfile(dateFolder, 'case_I_Te.png'), 'png');
% fprintf("Teのグラフを %s に保存しました。\n", fullfile(dateFolder, 'case_I_Te.png'));



% % neプロット作成
figure;
contourf(T_grid, R_grid, ne_data_2D, 1000, 'LineColor', 'none');
clim([0, 20e+19]); % カラーマップの範囲設定
xlim([475, 520]);
% カラーバーを追加してタイトルを設定
cb = colorbar;
cb.Title.String = 'n_{e} [m^{-3}]';
% 各要素の文字サイズを設定
xlabel('Time [\mus]', 'FontSize', 12);     % X軸ラベル
ylabel('R [m]', 'FontSize', 14); % Y軸ラベル
title('Case (a) plot of ne [m^{-3}]', 'FontSize', 16); % グラフタイトル
cb.Title.FontSize = 10;                    % カラーバーのタイトル文字サイズ
cb.FontSize = 10;                          % カラーバーの目盛りフォントサイズ
% 軸目盛りフォントサイズを設定
ax = gca; % 現在の座標軸を取得
ax.FontSize = 12; % 軸目盛りフォントサイズを設定
% カラーマップ設定
colormap('jet');
ylim([0.1, 0.3]);
xlim([470, 510]); 

% 保存する
% saveDir = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\figure';
% dateFolder = fullfile(saveDir, num2str(date));
% if ~exist(dateFolder, 'dir')
%     mkdir(dateFolder);
% end
% saveas(gcf, fullfile(dateFolder, 'case_I_ne.png'), 'png');
% fprintf("Teのグラフを %s に保存しました。\n", fullfile(dateFolder, 'case_I_ne.png'));

%%
% % ----------------------------
% % エラーバー付きプロット
% % ----------------------------
% selected_r = 2.5; % プロットしたいR値
% idx_selected = (r_list == selected_r);
% 
% % 平均値と標準偏差を計算
% Te_mean = mean(Te_matrix(idx_selected, :), 1, 'omitnan');
% Te_std = std(Te_matrix(idx_selected, :), 0, 1, 'omitnan');
% ne_mean = mean(ne_matrix(idx_selected, :), 0, 1, 'omitnan');
% ne_std = std(ne_matrix(idx_selected, :), 0, 1, 'omitnan');
% 
% % Teのエラーバー付きプロット
% figure;
% errorbar(time_us, Te_mean ./ K2ev, Te_std./ K2ev);
% xlim([450, 500])
% xlabel('Time [\mus]');
% ylabel('Te [eV]');
% title(['Te vs Time at R = ', num2str(selected_r)]);
% grid on;
% 
% % neのエラーバー付きプロット
% figure;
% errorbar(time_us, ne_mean, ne_std);
% xlabel('Time [\mus]');
% ylabel('ne [m^{-3}]');
% title(['ne vs Time at R = ', num2str(selected_r)]);
% grid on;
% 
% 
% ```