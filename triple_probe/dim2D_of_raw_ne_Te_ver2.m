%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Case-I, Case-Oのトリプルプローブデータから
% 2次元のTeとneを計算し、プロットする
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

% ファイルのパス設定
date = 240610; % case-I
% date = 240611; % case-O
filepath = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data\\';

% shotlist を修正 (1つ目の関数コードに合わせる)
% 1つ目の関数コードでは、dateが240610の時と240611の時でshotlistが切り替わっているので、それに合わせます。
if date == 240610
    shotlist = [18, 21, 24, 27, 33, 35, 38, 41]; % case-I
elseif date == 240611
    shotlist = [46, 49, 52, 55, 58, 61, 66, 69]; % case-O
else
    % デフォルトまたはエラー処理
    error('Unknown date selected for shotlist.');
end

% スプレッドシートからR方向データを取得
DOCID = '1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw'; % スプレッドシートのID
T = getTS6log(DOCID);
node = 'date';
T = searchlog(T, node, date);

% ショット番号に対応するインデックスを取得
shot_indices = arrayfun(@(x) find(T.shot == x, 1), shotlist);

% r_list を取得し、サイズを shotlist と一致させる
r_list = T.tripleProbeRPosition_cm_(shot_indices) * 1e-2; % R座標[m]

% パラメータ設定 (1つ目の関数に合わせる)
V2 = 20; % 電位設定
V3 = 10;
A = 1.0;  % Hガスの原子量
S_probe = 1.4e-5; % 電極表面積 [m^2]
current_noise_threshold = 1e-6; % 例: 1 microAmpere (A)以下はノイズと見なす
lowpass_freq = 0.1e6; % ローパスフィルターのカットオフ周波数
sampling_freq = 1e7; % サンプリング周波数
smoothing_window = 10; % 電流データのスムージングウィンドウ
final_smoothing_window = 7; % Te/ne行列の時間方向最終スムージングウィンドウ
% time_us_unifiedはここで定義
time_us = linspace(470, 550, 300); % 統一された時間軸 (470-550 μs, 300点)
% time_sliceはtime_rawの範囲を固定
time_slice = 4700:5500; % CSVファイルから読み込む行のインデックス

% 物理定数
q = 1.60217663e-19; % 電子電荷
kb = 1.380649e-23; % ボルツマン定数
K2ev = 11604.5250061657; % kelvin to eV
mi = A * 1.66054e-27; % イオン質量

% 時間データの統一と行列初期化
num_shots = length(shotlist);
Te_matrix = zeros(num_shots, length(time_us)); % Teデータ格納
ne_matrix = zeros(num_shots, length(time_us)); % neデータ格納
valid_shots_mask = false(1, num_shots); % valid_shots_mask を初期化

% シンボリック変数宣言
syms I1 I2 I3 Te ne_sym % ne はシンボリック変数なので ne_sym に変更

% Te解法
eqn_Te = (I1 + I2) / (I1 + I3) == ...
    (1 - exp(-q * V2 / (kb * Te))) / (1 - exp(-q * V3 / (kb * Te)));
Te_symbolic = solve(eqn_Te, Te);

% ne解法
eqn_ne = exp(-0.5) * S_probe * ne_sym * q * sqrt(kb * Te / mi) == ...
    (I3 - I2 * exp(-q * (V3 - V2) / (kb * Te))) / (1 - exp(-q * (V3 - V2) / (kb * Te)));
ne_symbolic = solve(eqn_ne, ne_sym);

% --- calculate_Te 関数の定義 (スクリプト内に含めるか、別ファイルにする) ---
% 今回はスクリプト内に含めます
function Te_values_K = calculate_Te(Te_symbolic, I1, I2, I3, K2ev_val)
    I1_c = I1; I2_c = I2; I3_c = I3;
    ratio_arg = (I1_c + I3_c) ./ (I2_c - I3_c);
    
    % problematic_indices_I2_minus_I3
    problem_idx = abs(I2_c - I3_c) < 1e-12 | ratio_arg <= 0 | abs(ratio_arg - 1) < 1e-12;
    I1_c(problem_idx) = NaN;
    I2_c(problem_idx) = NaN;
    I3_c(problem_idx) = NaN;
    
    Te_values_K = real(double(subs(Te_symbolic, {'I1', 'I2', 'I3'}, {I1_c, I2_c, I3_c})));
    Te_values_K(Te_values_K < 0 | Te_values_K > 100 * K2ev_val) = NaN; % K2evをK2ev_valに変更
end
% --- calculate_Te 関数定義終わり ---


% ショットごとのデータ読み込みと計算
for idx = 1:num_shots
    shot = shotlist(idx);
    filename = fullfile(filepath, num2str(date), ['ES_', num2str(date), sprintf('%03d', shot), '.csv']); % fullfileを使用
    
    if ~exist(filename, 'file')
        warning('File %s not found. Skipping this shot.', filename);
        continue; % ファイルが存在しない場合はこのショットをスキップ
    end
    
    valid_shots_mask(idx) = true; % 有効なショットをマーク
    
    % データ読み込み
    data = readmatrix(filename);

    time_raw = data(time_slice, 1); % 時間データ
    I2_raw = data(time_slice, 36); % 列36のデータ
    I3_raw = data(time_slice, 37); % 列37のデータ
    
    % Te/ne計算のためのデータ前処理
    % 負の値をNaNに置換
    I2_raw(I2_raw < 0) = NaN;
    I3_raw(I3_raw < 0) = NaN;
    
    % NaNを線形補間し、lowpassフィルタを適用
    I2_processed = lowpass(fillmissing(I2_raw, 'linear'), lowpass_freq, sampling_freq);
    I3_processed = lowpass(fillmissing(I3_raw, 'linear'),lowpass_freq, sampling_freq);
    
    % フィルタ後のデータに、元のNaNの位置を復元
    I2_processed(isnan(I2_raw)) = NaN;
    I3_processed(isnan(I3_raw)) = NaN;

    % 信号が小さすぎる区間を無効化するための閾値処理
    peak_current = max(max(I2_processed, [], 'omitnan'), max(I3_processed, [], 'omitnan'));
    if ~isempty(peak_current) && peak_current > 0
        current_threshold = peak_current * 0.05; % ピーク値の5%を閾値とする
        I2_processed(I2_processed < current_threshold) = NaN;
        I3_processed(I3_processed < current_threshold) = NaN;
    end
    
    I1_processed = I2_processed + I3_processed;

    % Te 計算 (calculate_Te 関数を使用)
    try
        Te_values_K = calculate_Te(Te_symbolic, I1_processed, I2_processed, I3_processed, K2ev); % K2evを渡す
    catch ME
        disp(['An error occurred during Te calculation for shot ', num2str(shot), ':']);
        disp(ME.message);
        Te_values_K = NaN(size(I1_processed));
    end

    % ne 計算
    try
        ne_values = real(double(subs(ne_symbolic, {Te, I2, I3}, {Te_values_K, I2_processed, I3_processed})));
        ne_values(ne_values < 0) = NaN;
    catch ME
        disp(['An error occurred during ne calculation for shot ', num2str(shot), ':']);
        disp(ME.message);
        ne_values = NaN(size(I2_processed));
    end

    % 補間して時間軸を統一 (NaNをNaNで外挿, 1つ目のコードのcalculate_Teの後にNaNを埋めてないのでこのままで良い)
    Te_matrix(idx, :) = interp1(time_raw, Te_values_K, time_us, 'linear', NaN); % extrapval を NaN に変更
    ne_matrix(idx, :) = interp1(time_raw, ne_values, time_us, 'linear', NaN); % extrapval を NaN に変更

    % 各ショットのTe/neデータ（時間軸方向）の最終スムージング
    % % これはTe_matrix_shotsに格納される前に適用される (1つ目の関数に合わせる)
    % Te_matrix(idx, :) = smoothdata(Te_matrix(idx, :), 'movmedian', final_sm_win, 'omitnan');
    % ne_matrix(idx, :) = smoothdata(ne_matrix(idx, :), 'movmedian', final_sm_win, 'omitnan');

    % plot_shotlist の定義
    plot_shotlist = []; % プロットしたいショット番号を直接指定
    if ismember(shot, plot_shotlist)
        figure('Name', ['Raw I2 and I3 for Shot ', num2str(shot)]);
        hold on;
        grid on;
        % 生信号
        plot(time_raw, I2_raw, 'b-', 'DisplayName', 'I2 (Raw)'); % I2_rawを使用
        plot(time_raw, I3_raw, 'r-', 'DisplayName', 'I3 (Raw)'); % I3_rawを使用
        % 補間後のデータ (Te, ne計算に使用された範囲)
        plot(time_raw, I2_processed, 'c-', 'LineWidth', 1.5, 'DisplayName', 'I2 (Corrected)'); % I2_processedを使用
        plot(time_raw, I3_processed, 'm-', 'LineWidth', 1.5, 'DisplayName', 'I3 (Corrected)'); % I3_processedを使用
        % 閾値
        yline(current_threshold, 'k--', 'LineWidth', 1, 'DisplayName', 'Noise Threshold');
        xlim([450, 520]);
        xlabel('Time [μs]');
        ylabel('Current [A]');
        title(['Raw I2 and I3 for Shot ', num2str(shot), ' (R = ', num2str(r_list(idx)), ' m)']);
        legend show;
        hold off;

        % CurrentRatio のプロット
        CurrentRatio = (I1_processed + I2_processed) ./ (I1_processed + I3_processed); % Te式に基づく比率の分子はI1+I2, 分母はI1+I3。電流比だけならI1+I2 / I1+I3 でいいが、電流値の比ならI2/I3などが考えられる。ここではTe計算に使ったI1,I2,I3_processedを使用
        figure; plot(time_raw, CurrentRatio); title(['Current Ratio for Shot ', num2str(shot)]);
        xlim([450, 520]);
        ylim([-5, 5]); % 適宜範囲を調整
        grid on;
    end
end

%%
% --- 集計処理 ---
% 有効なショットのデータのみを抽出 (valid_idxではなくvalid_shots_maskを使用)
Te_matrix = Te_matrix(valid_shots_mask, :);
ne_matrix = ne_matrix(valid_shots_mask, :);
r_list_valid = r_list(valid_shots_mask);

% 時間方向のNaNが残る場合はここでfillmissing (Te_matrix_shotsに相当する部分)
% （Te_matrix内のsmoothdataでomitnanを使っているため、通常はNaNが残らないはずですが、念のため）
% for i = 1:size(Te_matrix, 1)
%     Te_matrix(i, :) = fillmissing(Te_matrix(i, :), 'linear', 'EndValues', 'nearest'); % 端のNaNも埋める
%     ne_matrix(i, :) = fillmissing(ne_matrix(i, :), 'linear', 'EndValues', 'nearest');
% end

% 各ショットの時間データについて、内側のNaNのみを線形補間します。
% 'EndValues'オプションを削除することで、始点・終点のNaNはそのまま残します。
for i = 1:size(Te_matrix, 1)
    Te_matrix(i, :) = fillmissing(Te_matrix(i, :), 'linear', 'EndValues', 'none');
    ne_matrix(i, :) = fillmissing(ne_matrix(i, :), 'linear', 'EndValues', 'none');
end


% 同じr_listの値でTeとneの平均値を計算
unique_r = unique(r_list_valid); % r_list_validを使用
Te_avg_matrix = zeros(length(unique_r), length(time_us));
ne_avg_matrix = zeros(length(unique_r), length(time_us));

for i = 1:length(unique_r)
    idx_r = (r_list_valid == unique_r(i)); % r_list_validを使用
    Te_avg_matrix(i, :) = mean(Te_matrix(idx_r, :), 1, 'omitnan'); % NaN は 'omitnan' で無視
    ne_avg_matrix(i, :) = mean(ne_matrix(idx_r, :), 1, 'omitnan'); % NaN は 'omitnan' で無視
end

% 平均化後の時間方向の最終スムージング
for i = 1:size(Te_avg_matrix, 1)
    Te_avg_matrix(i, :) = smoothdata(Te_avg_matrix(i, :), 'movmedian', 7, 'omitnan'); 
    ne_avg_matrix(i, :) = smoothdata(ne_avg_matrix(i, :), 'movmedian', 7, 'omitnan');
end

% 各時刻のデータ点について、R方向（列方向）のNaNを線形補間します。
% これにより、特定の半径でデータが全くなかった場合も、その上下のデータから補間されます。
Te_avg_matrix = fillmissing(Te_avg_matrix, 'linear', 1, 'EndValues', 'none');
ne_avg_matrix = fillmissing(ne_avg_matrix, 'linear', 1, 'EndValues', 'none');

% さらに、R方向に平滑化をかけ、より滑らかな分布にします。
% ウィンドウサイズは時間方向より小さい値(例: 3)が適切です。
Te_avg_matrix = smoothdata(Te_avg_matrix, 1, 'movmedian', 3, 'omitnan');
ne_avg_matrix = smoothdata(ne_avg_matrix, 1, 'movmedian', 3, 'omitnan');

% 平均化後の最終スムージング (1つ目の関数に合わせる)
% ここでスムージングをかけない場合はこれらの行をコメントアウト
% for i = 1:size(Te_avg_matrix, 1) % 時間方向
%     Te_avg_matrix(i, :) = smoothdata(Te_avg_matrix(i, :), 'movmedian', 7, 'omitnan'); 
%     ne_avg_matrix(i, :) = smoothdata(ne_avg_matrix(i, :), 'movmedian', 7, 'omitnan');
% end


% --------------------------------------------------------------------------------------
% ここから先はMATファイル作製とコンタープロット（ほぼ変更なし）
% --------------------------------------------------------------------------------------

% 必要なデータを準備
R_values = unique_r;             % R方向の値
time_values = time_us;           % 時間軸
Te_data_2D = Te_avg_matrix ./ K2ev; % Teデータを[eV]に変換
ne_data_2D = ne_avg_matrix;       % neデータ

% 3次元データ構造を作成
triple_data2D = struct(); 
triple_data2D.Te_avg = reshape(Te_data_2D, [length(R_values), 1, length(time_values)]);
triple_data2D.ne_avg = reshape(ne_data_2D, [length(R_values), 1, length(time_values)]);


% MATファイルとして保存 (1つ目の関数に合わせるために一部変更)
targetFolderPath = fullfile(filepath, '..', 'triple_probe', 'mat_test', num2str(date)); % 相対パス調整

if ~exist(targetFolderPath, 'dir')
    mkdir(targetFolderPath);
end

% ファイル名の動的決定 (1つ目の関数のcase_name引数に相当するものをスクリプトで設定)
if date == 240610
    case_name = 'case-I';
elseif date == 240611
    case_name = 'case-O';
else
    case_name = 'unknown_case';
end
matFileName = fullfile(targetFolderPath, ['triple_data2D_', case_name, '.mat']);

% 保存する変数リストも1つ目の関数に合わせる
save(matFileName, ...
    'triple_data2D', ...
    'R_values', ...
    'time_values', ...
    'shotlist', ... % 元のshotlistを保存
    'r_list_valid', ... % valid_shots_maskを適用したr_listを保存
    'time_raw', ...     % 各ショット共通のtime_rawを保存 (もし毎回同じなら)
    'I2_raw', ...       % 各ショットの最後のI2_raw_all_shots (または適切なもの)
    'I3_raw', ...       % 各ショットの最後のI3_raw_all_shots (または適切なもの)
    'time_slice');
disp(['MATファイル ', matFileName, ' を保存しました。']);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% conturfのプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATファイルを読み込む
load(matFileName); 

% 必要なデータを取り出す
Te_data_2D = squeeze(triple_data2D.Te_avg(:, 1, :)); % 'Te' から 'Te_avg' に変更
ne_data_2D = squeeze(triple_data2D.ne_avg(:, 1, :)); % 'ne' から 'ne_avg' に変更

% グリッドデータを作成
[T_grid, R_grid] = meshgrid(time_values, R_values); 

% Teプロット作成
figure;
contourf(T_grid, R_grid, Te_data_2D, 10000, 'LineColor', 'none');
clim([0, 15]); 
xlim([470, 510]); % 修正
ylim([0.1, 0.3]); % 修正
% カラーバーを追加してタイトルを設定
cb = colorbar;
cb.Title.String = 'T_{e} [eV]';
% 個別に文字サイズを設定
xlabel('Time [\mus]', 'FontSize', 12);     
ylabel('R [m]', 'FontSize', 14);           
title(['Case - ', case_name, ' plot of Te [eV]'], 'FontSize', 16); % タイトル修正
cb.Title.FontSize = 10;                    
cb.FontSize = 10;                          
% 軸目盛りフォントサイズを設定
ax = gca; 
ax.FontSize = 12; 
% カラーマップ設定
colormap('jet');

% neプロット作成
figure;
contourf(T_grid, R_grid, ne_data_2D, 100, 'LineColor', 'none');
clim([0, 1.5e+20]); 
xlim([470, 510]); % 修正
ylim([0.1, 0.3]); % 修正
% カラーバーを追加してタイトルを設定
cb = colorbar;
cb.Title.String = 'n_{e} [m^{-3}]';
% 各要素の文字サイズを設定
xlabel('Time [\mus]', 'FontSize', 12);     
ylabel('R [m]', 'FontSize', 14); 
title(['Case - ', case_name, ' plot of ne [m^{-3}]'], 'FontSize', 16); % タイトル修正
cb.Title.FontSize = 10;                    
cb.FontSize = 10;                          
% 軸目盛りフォントサイズを設定
ax = gca; 
ax.FontSize = 12; 
% カラーマップ設定
colormap('jet');