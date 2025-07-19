%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2次元プロットして、確認する
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% クリアとデータの読み込み
clear all;

% データを読み込む
% load("G:\My Drive\lab\lab_data\mach_probe_rawdata\240611\070.mat");
load("G:\My Drive\lab\lab_data\mach_probe_rawdata\240610\047_4700-5200.mat");

% 水素イオンの質量 (mi) を設定 [kg]
mi = 1.6726219e-27; % 水素イオンの質量

% 2次元カウンタープロットの作成
if isfield(MPdata2D, 'flow_forplot') && isfield(MPdata2D, 'zq') && isfield(MPdata2D, 'rq')
    % flow_forplot データを取り出し
    v = squeeze(MPdata2D.flow_forplot(151, :, :)); % flow_forplot の速度成分
    v = abs(v);

    % カウンタープロット
    figure;
    contourf(MPdata2D.zq, MPdata2D.rq, v, 201, 'edgecolor', 'none'); % 100等高線レベル
    colorbar;
    xlabel('Z-axis [m]');
    ylabel('R-axis [m]');
    title('2D Contour Plot of Velocity (flow\_forplot)');
else
    disp('MPdata2D does not contain necessary fields for contour plot.');
end

% 運動エネルギーの計算
if exist('v', 'var')
    % 運動エネルギーの計算 (U = 1/2 * mi * v^2)
    U = 0.5 * mi .* (v.^2); % 要素ごとの計算

    % NaN を無視して運動エネルギーの総和を計算
    U(isnan(U)) = 0; % NaN を 0 に置き換える
    % 運動エネルギーの総和を計算
    total_kinetic_energy = sum(U(:)); % すべての値を足し算

    % 結果を表示
    fprintf('運動エネルギーの総和: %.3e J\n', total_kinetic_energy);

else
    disp('Velocity data (flow_forplot) is not available. Unable to compute kinetic energy.');
end


%%
% 水素イオンの質量 (mi) を設定 [kg]
mi = 1.6726219e-27; % 水素イオンの質量
ni = 1e20;

% 運動エネルギーを格納する配列
kinetic_energy_array = [];

% trangeの範囲を設定 (例: 1から100)
user_defined_trange = 1:501; % 必要に応じて変更してください

% データが正しい形式で存在するか確認
if isfield(MPdata2D, 'flow_forplot') && isfield(MPdata2D, 'trange') ...
        && isfield(MPdata2D, 'zq') && isfield(MPdata2D, 'rq')

    % ユーザー定義のtrangeに基づいて運動エネルギーを計算
    for t_idx = user_defined_trange
        % t_idxがデータ範囲を超えていないか確認
        if t_idx > length(MPdata2D.trange)
            disp('Warning: 指定した範囲がデータの範囲を超えています。計算をスキップします。');
            continue;
        end

        % 速度データ (flow_forplot) を取り出し
        v = squeeze(MPdata2D.flow_forplot(t_idx, :, :)); % 現在の時間フレームの速度データ
        v = abs(v) * 1e3; % 絶対値を取り、km/s を m/s に変換

        % 運動エネルギーの計算 (U = 1/2 * mi * v^2)
        U = 0.5 * mi .* ni .* (v.^2); % 要素ごとの計算

        % NaN を無視して総和を計算
        U(isnan(U)) = 0; % NaN を 0 に置き換える
        total_energy_t = sum(U(:)); % 現在の時間フレームの運動エネルギーを計算

        % 配列に追加
        kinetic_energy_array = [kinetic_energy_array, total_energy_t];
    end

    % 全フレームの運動エネルギーの総和を計算
    total_kinetic_energy = sum(kinetic_energy_array);

    % 結果を表示
    fprintf('各時間フレームの運動エネルギーの総和: %.3e J\n', total_kinetic_energy);

else
    disp('MPdata2D does not contain necessary fields for computation.');
end

%%
plot(user_defined_trange, kinetic_energy_array)
scatter(user_defined_trange, kinetic_energy_array)