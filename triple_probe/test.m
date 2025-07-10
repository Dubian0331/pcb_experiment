% 定数を定義
q = 1.6e-19; % 電子の電荷 (クーロン)
kB = 1.38e-23; % ボルツマン定数 (ジュール/ケルビン)

% 変数を数値で定義
I1_val = 1e-3; % アンペア
I2_val = 2e-3; % アンペア
I3_val = 3e-3; % アンペア
V2_val = 10;    % ボルト
V3_val = 20;   % ボルト

% 最適化関数を定義
objective = @(Te) abs(((I1_val + I2_val) / (I1_val + I3_val)) - ...
    ((1 - exp(-q * V2_val / (kB * Te))) / (1 - exp(-q * V3_val / (kB * Te)))));

% 初期値と範囲を設定
Te_initial = 1; % 初期値 (ケルビン)
Te_lower_bound = 1e-2; % 下限値 (ケルビン)
Te_upper_bound = 1e4; % 上限値 (ケルビン)

% 最適化を実行
options = optimset('Display', 'iter', 'TolX', 1e-6);
[Te_optimal_K, fval] = fminbnd(objective, Te_lower_bound, Te_upper_bound, options);

% K から eV への変換
Te_optimal_eV = Te_optimal_K * kB / q;

% 結果を表示
fprintf('Teの最適解は %.5f eV です。\n', Te_optimal_eV);

% グラフを作成
Te_values_K = linspace(Te_lower_bound, Te_upper_bound, 1000);
Te_values_eV = Te_values_K * kB / q;
objective_values = arrayfun(objective, Te_values_K);

figure;
plot(Te_values_eV, objective_values, 'LineWidth', 2);
hold on;
plot(Te_optimal_eV, fval, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

title('Objective Function vs Te');
xlabel('Te (eV)');
ylabel('Objective Function Value');
grid on;


%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1次元のプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% データの読み込み
date = 240610;
prompt = 'Enter the shot number: ';
dlgtitle = 'Shot Number';
dims = [1 30];
definput = {''};
shotnum = inputdlg(prompt, dlgtitle, dims, definput);
shotnum = str2double(shotnum{1});

if shotnum < 10
    shotnum = ['00', num2str(shotnum)];
else 
    shotnum = ['0', num2str(shotnum)];
end

save_filepath = "G:\My Drive\lab\lab_data\mach_probe_rawdata";

% 研究室内にいる場合のファイルパス
% filepath = "\\NIFS\experiment\results\MachProbe\";
% 研究室外のファイルパス
filepath = "G:\My Drive\lab\lab_data\mach_probe_rawdata";

filename = strcat(filepath, '\', num2str(date), '\ES_', num2str(date), shotnum, '.csv');

% ファイルの読み込み
test = readmatrix(filename);
index_start = 4500;
index_end = 5000;
V2 = 20;
V3 = 10;

% if you use H gas. The A is 1.00.
A = 1.00; % atomic weight
% if you use Ar gas. The A is 39.95.
% A = 39.95;
S_probe = 2.89*10^-5; % probe surface area [m2] Φ0.45mm electrode
I2_values = test(index_start:index_end, 36);
I3_values = test(index_start:index_end, 37);

% 強化したスムージング
I2_values = smoothdata(I2_values, 'movmean', 10);
I3_values = smoothdata(I3_values, 'movmean', 10);

I1_values = I2_values + I3_values;
time = test(index_start:index_end, 1);

q = 1.60217663 * 10^(-19); % electron charge
kB = 1.380649 * 10^(-23); % boltzmann const
K2ev = 11604.5250061657; % kelvin to eV
mi = A * 1.66054 * 10^(-27); % ion mass

% solve for Te using optimization
Te_values_eV = zeros(size(I1_values));
options = optimset('Display', 'off', 'TolX', 1e-6);
for i = 1:length(I1_values)
    I1 = I1_values(i);
    I2 = I2_values(i);
    I3 = I3_values(i);
    if I1 > 0 && I3 > 0
        objective = @(Te) abs(((I1 + I2) / (I1 + I3)) - ...
            ((1 - exp(-q * V2 / (kB * Te))) / (1 - exp(-q * V3 / (kB * Te)))));
        Te_optimal_K = fminbnd(objective, 0.1, 10e+10, options); % Set a wider range
        Te_values_eV(i) = Te_optimal_K * kB / q;
    else
        Te_values_eV(i) = NaN; % Handle invalid data points
    end
end

% solve for ne
syms ne Te I2 I3
eqn = exp(-0.5) * S_probe * ne * q * (kB * Te / mi)^0.5 == (I3 - I2*exp(-q*(V3-V2)/(kB*Te)))/(1-exp(-q*(V3-V2)/(kB*Te)));
ne_symbolic = solve(eqn,ne);
% Calculate ne_values and check for division by zero
try
    ne_values = real(double(subs(ne_symbolic, {Te I2 I3}, {Te_values_eV I2_values I3_values})));
catch
    % Handle division by zero error
    warning('Division by zero occurred while solving for ne. Applying previous valid value as fallback.');
    ne_values = zeros(size(I2_values));
    for i = 2:length(I2_values)
        try
            ne_values(i) = real(double(subs(ne_symbolic, {Te I2 I3}, {Te_values_eV(i) I2_values(i) I3_values(i)})));
        catch
            if i > 1
                ne_values(i) = ne_values(i - 1);
            end
        end
    end
end

% plot
figure;
sgtitle('Triple probe data of H spheromak')

subplot(3,1,1)
hold on;
plot(time, I2_values);
plot(time, I3_values);
legend({['I2 (V_p=', num2str(V2), 'V)'], ['I3 (V_p=', num2str(V3), 'V)']})
ylabel('Probe current [A]')
xlim([450, 500])
hold off
grid on

subplot(3,1,2)
scatter(time, Te_values_eV, 5, "black");
xlim([450, 500]);
ylim([0 20]);
ylabel('Te [eV]')
grid on

subplot(3,1,3)
scatter(time, ne_values, 3, "black");
xlim([450, 500]);
xlabel('time [us]')
ylabel('ne [m^{-3}]')
grid on

mkdir(strcat(save_filepath, '\', num2str(date), '\figure\triple_probe'));
saveas(gcf, strcat(save_filepath, '\', num2str(date), '\figure\triple_probe\', shotnum, '_', num2str(index_start), '-', num2str(index_end), 'us'), 'png');
fprintf("save your file %s", shotnum);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2次元プロット - グリッド化による最適化適用版
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ファイルのパス設定
date = 240610; % 日付
filepath = "G:\My Drive\lab\lab_data\mach_probe_rawdata\"; % ファイルパス
shotlist = 17:47; % 読み込むショット番号リスト

% スプレッドシートからR方向データを取得
DOCID = '1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw'; % スプレッドシートのID
T = getTS6log(DOCID);
node = 'date';
T = searchlog(T, node, date);

% ショット番号に対応するインデックスを取得
shot_indices = arrayfun(@(x) find(T.shot == x, 1), shotlist);

% r_list を取得し、サイズを shotlist と一致させる
r_list = T.tripleProbeRPosition_cm_(shot_indices) * 1e-1; % R座標[m]

% パラメータ設定
V2 = 20;
V3 = 10; % 電位設定
A = 1.0; % Hガスの場合の原子量
S_probe = 2.89e-5; % 電極表面積 [m^2]
q = 1.60217663e-19; % 電子電荷
kb = 1.380649e-23; % ボルツマン定数
mi = A * 1.66054e-27; % イオン質量
K2ev = 11604.5250061657; % kelvin to eV

% 時間データの統一と行列初期化
num_shots = length(shotlist);
time_us = linspace(450, 600, 300); % 統一された時間軸 (450-600 µs, 300点)
Te_matrix = zeros(num_shots, length(time_us)); % Teデータ格納
ne_matrix = zeros(num_shots, length(time_us)); % neデータ格納

% Te グリッドの設定
Te_grid_K = linspace(0.1, 30000, 5000); % Kelvin グリッド
Te_grid_eV = Te_grid_K * kb / q; % eV グリッド

% シンボリック変数宣言
syms ne Te I2 I3

% ne の方程式
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
        time = data(4500:7500, 1); % 時間データ
        I2_values = data(4500:7500, 36); % 列36のデータ
        I3_values = data(4500:7500, 37); % 列37のデータ
        
        % 欠損データの補完とスムージング
        I2_values = fillmissing(I2_values, 'linear');
        I3_values = fillmissing(I3_values, 'linear');
        I2_values = smoothdata(I2_values, 'movmean', 10);
        I3_values = smoothdata(I3_values, 'movmean', 10);
        I1_values = I2_values + I3_values;

        % Te のグリッド計算
        Te_values_eV = zeros(size(I1_values));
        for i = 1:length(I1_values)
            I1 = I1_values(i);
            I2 = I2_values(i);
            I3 = I3_values(i);
            if I1 > 0 && I3 > 0
                % グリッド上で目的関数を評価
                objective_values = abs(((I1 + I2) / (I1 + I3)) - ...
                    ((1 - exp(-q * V2 ./ (kb * Te_grid_K))) ./ ...
                    (1 - exp(-q * V3 ./ (kb * Te_grid_K)))));
                [~, min_idx] = min(objective_values);
                Te_values_eV(i) = Te_grid_eV(min_idx);
            else
                Te_values_eV(i) = NaN;
            end
        end

        % ne 計算
        ne_values = zeros(size(Te_values_eV));
        for i = 1:length(Te_values_eV)
            try
                if ~isnan(Te_values_eV(i)) && Te_values_eV(i) > 0
                    ne_values(i) = real(double(subs(ne_symbolic, ...
                        {Te, I2, I3}, {Te_values_eV(i), I2_values(i), I3_values(i)})));
                else
                    ne_values(i) = NaN;
                end
            catch
                ne_values(i) = NaN;
            end
        end

        % 補間して時間軸を統一
        Te_matrix(valid_idx, :) = interp1(time, Te_values_eV, time_us, 'linear', NaN);
        ne_matrix(valid_idx, :) = interp1(time, ne_values, time_us, 'linear', NaN);
    else
        warning('File %s not found. Skipping this shot.', filename);
    end
end

% 有効データのみを使用
Te_matrix = Te_matrix(1:valid_idx, :);
ne_matrix = ne_matrix(1:valid_idx, :);
r_list = r_list(1:valid_idx);

% ----------------------------
% meshgrid を用いた補間とプロット
% ----------------------------
unique_r = unique(r_list); % 一意なR値
[X, Y] = meshgrid(time_us, unique_r); % 元のグリッド
[Xq, Yq] = meshgrid(linspace(min(time_us), max(time_us), 100), ...
                    linspace(min(unique_r), max(unique_r), 100)); % 補間用グリッド

% Te の補間
Te_smooth = griddata(repmat(time_us, valid_idx, 1), repmat(r_list, 1, length(time_us)), Te_matrix, Xq, Yq, 'cubic');
figure;
imagesc(Xq(1, :), Yq(:, 1), Te_smooth);
set(gca, 'YDir', 'normal');
colorbar;
clim([0, 15]); % Te の表示範囲 (eV)
xlabel('Time [\mus]');
ylabel('R position [m]');
title('Smoothed 2D Plot of Te [eV]');
colormap(jet);

% ne の補間
ne_smooth = griddata(repmat(time_us, valid_idx, 1), repmat(r_list, 1, length(time_us)), ne_matrix, Xq, Yq, 'cubic');
figure;
imagesc(Xq(1, :), Yq(:, 1), ne_smooth);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('Time [\mus]');
ylabel('R position [m]');
title('Smoothed 2D Plot of ne [m^{-3}]');
colormap(jet);
