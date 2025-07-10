function plot_iv_curve_and_calculate_te(date, shotnum, V)
    % 引数:
    %   date: 日付（例: 240608）
    %   shotnum: ショット番号のリスト（例: [2, 5, 8, 9, 10]）
    %   V: バイアス電圧の配列
    
    % 指定された時間のリスト
    times_us = [466, 470, 480, 490, 500];
    num_times = length(times_us);
    
    % カラーマップの設定
    colors = lines(num_times);
    markers = {'o', 'x', 's', 'd', '^'};
    
    % プロットの初期化
    figure;
    hold on;
    
    % データ格納用配列の初期化
    I_all = zeros(length(V), num_times);
    
    % 各ショット番号について処理を実行
    for i = 1:length(shotnum)
        % ショット番号を適切な形式に変換
        if shotnum(i) < 10
            shot_str = strcat('00', num2str(shotnum(i)));
        else
            shot_str = strcat('0', num2str(shotnum(i)));
        end
        
        % ファイルパスの生成
        file_path = strcat('\\NIFS\\experiment\\results\\MachProbe\\', num2str(date), '\\ES_', num2str(date), shot_str, '.csv');
        
        % CSVファイルからデータを読み込む
        data = readmatrix(file_path);
        
        % 時間の列を抽出
        time_column = data(:, 1); % 最初の列が時間データと仮定
        
        % 各指定時間の電流データを抽出（例としてch37列を使用）
        for j = 1:num_times
            % 指定時間に最も近いデータポイントを見つける
            [~, idx] = min(abs(time_column - times_us(j)));
            I_all(i, j) = data(idx, 37); % 37列目がch9-4のデータと仮定
        end
    end
    
    % データをプロット
    for j = 1:num_times
        plot(V, I_all(:, j), strcat('-', markers{j}), 'Color', colors(j, :), 'LineWidth', 2, 'MarkerSize', 6, 'DisplayName', sprintf('%dus', times_us(j)));
    end
    
    % 微分計算
    dIdV = zeros(size(I_all));
    for j = 1:num_times
        dIdV(:, j) = gradient(I_all(:, j), V);
    end
    
    % 微分値の取得
    dIdV_0 = mean(dIdV(1:2, :)); % V=0 付近の微分値
    dIdV_inf = mean(dIdV(end-2:end, :)); % V>>0 付近の微分値

    % I0 の計算 (dIdV_0 と dIdV_inf を用いて接線の交点を求める)
    V0 = 0; % V=0 付近での電圧
    V_inf = V(end); % V が大きいときの電圧（ここでは最大電圧を使用）
    I0 = mean(I_all(1:2, :)) - dIdV_0 * V(1); % V=0 での電流値
    I_inf = mean(I_all(end-2:end, :)) - dIdV_inf * V(end); % V=最大電圧での電流値
    
    % 電子温度の計算
    k_e = 1.38e-23; % ボルツマン定数
    e = 1.602e-19; % 電子の電荷
    
    Te = - (I0 ./ (2 * dIdV_0 - 1.64 * dIdV_inf)) * (k_e / e);
    
    % 結果の表示
    fprintf('電子温度 (Te) の計算結果:\n');
    for j = 1:num_times
        fprintf('%dus: Te = %.2f eV\n', times_us(j), Te(j));
    end

    % 接線のプロット
    for j = 1:num_times
        % 接線のプロット用の電圧範囲
        V_line = linspace(min(V), max(V), 100);
        
        % dIdV_0 による接線
        I_line_0 = dIdV_0(j) * (V_line - V(1)) + mean(I_all(1:2, j));
        plot(V_line, I_line_0, '--', 'Color', colors(j, :));
        
        % dIdV_inf による接線
        I_line_inf = dIdV_inf(j) * (V_line - V(end)) + mean(I_all(end-2:end, j));
        plot(V_line, I_line_inf, ':', 'Color', colors(j, :));
        
        % I0 のプロット
        plot(0, I0(j), 'p', 'Color', colors(j, :), 'MarkerSize', 10, 'MarkerFaceColor', colors(j, :), 'DisplayName', sprintf('I0 %dus', times_us(j)));
    end
    
    % グラフの設定
    xlabel('bias voltage [V]');
    ylabel('probe current [A]');
    title('IV characteristic curve');
    legend show;
    grid on;
    hold off;
end

%%
date = 240609;
shotnum = [24, 25, 26, 35, 27, 28, 30, 36, 37, 38, 39, 40, 41];
V = [-40, -30, -20, -15, -10, -5, 0, 5, 10, 15, 20, 30, 40]; % バイアス電圧

plot_iv_curve_and_calculate_te(date, shotnum, V);






%%
date = 240609;
shotnum = [24, 25, 26, 35, 27, 28, 30, 36, 37, 38, 39, 40, 41];
V = [-40, -30, -20, -15, -10, -5, 0, 5, 10, 15, 20, 30, 40]; % バイアス電圧

plot_iv_curve_and_calculate_te(date, shotnum, V);
