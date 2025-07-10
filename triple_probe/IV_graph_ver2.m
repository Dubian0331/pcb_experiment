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
            I_all(i, j) = -data(idx, 37); % 37列目がch9-4のデータと仮定し、正負を逆にする
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
    
    % I0 の計算 (I0 は V が大きいときの電流の平均値と仮定)
    I0 = mean(I_all(end-3:end, :)); % V が大きいときの値を平均
    
    % 微分値の取得
    dIdV_0 = mean(dIdV(1:2, :)); % V=0 付近の微分値
    dIdV_inf = mean(dIdV(end-3:end, :)); % V>>0 付近の微分値
    
    % 電子温度の計算
    k_e = 1.38e-23; % ボルツマン定数
    e = 1.602e-19; % 電子の電荷
    
    Te = - (I0 ./ (2 * dIdV_0 - 1.64 * dIdV_inf)) * (k_e / e);
    
    % 結果の表示
    fprintf('電子温度 (Te) の計算結果:\n');
    for j = 1:num_times
        fprintf('%dus: Te = %.2f eV\n', times_us(j), Te(j));
    end
    
    % グラフの設定
    xlabel('bias voltage [V]');
    ylabel('probe current [A]');
    title('IV characteristic curve');
    legend('show', 'Location', 'northeastoutside');
    text(-35, max(max(I_all)) * 0.9, sprintf('Times: %dus', times_us), 'FontSize', 10, 'VerticalAlignment', 'top');
    grid on;
    hold off;
end

%%
date = 240609;
shotnum = [24, 25, 26, 35, 27, 28, 30, 36, 37, 38, 39, 40, 41];
V = [-40, -30, -20, -15, -10, -5, 0, 5, 10, 15, 20, 30, 40]; % バイアス電圧

plot_iv_curve_and_calculate_te(date, shotnum, V);
