%% ========================================================================
%  I_ratio の計算・保存・プロット用スクリプト (直接入力版)
% =========================================================================
clear;
% close all;

% パラメータ
date = 240611;
shotnum = 47;
K = 2.5;  % 比例定数
k = 1.38e-23;   % [J/K]
q = 1.60e-19;   % [C]
Te = 10 * q / k; % [K]　% T_e：電子温度→青山さんのp.35を参照 5~8[eV]
mi = 1.672e-27; % [kg]
cs = sqrt(2 * k * Te / mi); % [m/s]


shotnum_orig = shotnum; % プロットタイトル用に元のショット番号を保持

if shotnum < 10
    shotnum_str = ['00',num2str(shotnum)];
elseif shotnum < 100
    shotnum_str = ['0',num2str(shotnum)];
else 
    shotnum_str = num2str(shotnum);
end

% ファイルパスの構築
filepath = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\MachProbe_data';
filename = strcat(filepath,'\',num2str(date),'\ES_',num2str(date),shotnum_str,'.csv');
fprintf('読み込みファイル: %s\n', filename);

try
    df = readtable(filename, 'Delimiter', ',',"VariableNamingRule","preserve");
catch ME
    errordlg(sprintf('ファイルが読み込めませんでした。\nエラー: %s', ME.message),'ファイルエラー');
    return;
end

% I_ratioの計算
datalen = 10001;
time_vector = df.("time(us)")(1:datalen);

I_ratio_list = table(time_vector, 'VariableNames', {'time'});
Iup_list = table(time_vector, 'VariableNames', {'time'});
Idown_list = table(time_vector, 'VariableNames', {'time'});
flow_list = table(time_vector, 'VariableNames', {'time'});
Iup_raw_list = table(time_vector, 'VariableNames', {'time'});
Idown_raw_list = table(time_vector, 'VariableNames', {'time'});
current_threshold_list = zeros(1, 9);

channel_list = [1, 3];

% --- 各プローブヘッドのデータをループ処理 ---
for i = 1:9
    ch_up_name = sprintf('ch%d', 4*(i-1) + channel_list(1));
    ch_down_name = sprintf('ch%d', 4*(i-1) + channel_list(2));
    
    Iup_raw = df.(ch_up_name);
    Idown_raw = df.(ch_down_name);
    Iup_raw_list.(['ch' num2str(i)]) = Iup_raw(1:datalen);
    Idown_raw_list.(['ch' num2str(i)]) = Idown_raw(1:datalen);
    
    Iup_raw(Iup_raw < 0) = NaN;
    Idown_raw(Idown_raw < 0) = NaN;
    
    Iup_for_filter = fillmissing(Iup_raw, 'linear');
    Idown_for_filter = fillmissing(Idown_raw, 'linear');
    
    Iup_filtered = lowpass(Iup_for_filter, 1e5, 10e6);
    Idown_filtered = lowpass(Idown_for_filter, 1e5, 10e6);

    % 元のNaNの位置を復元し、不適切なデータ点を除外する
    Iup_processed = Iup_filtered;
    Idown_processed = Idown_filtered;
    Iup_processed(isnan(Iup_raw)) = NaN;
    Idown_processed(isnan(Idown_raw)) = NaN;

    % 処理後のIup, Idownを保存
    Iup_list.(['ch' num2str(i)]) = Iup_processed(1:datalen);
    Idown_list.(['ch' num2str(i)]) = Idown_processed(1:datalen);

    % 閾値を動的に決定
    peak_current = max(max(Iup_processed, [], 'omitnan'), max(Idown_processed, [], 'omitnan'));
    current_threshold = peak_current * 0.01;
    current_threshold_list(i) = current_threshold;

    % 閾値を下回るデータをNaNにする
    Iup_processed(Iup_processed < current_threshold) = NaN;
    Idown_processed(Idown_processed < current_threshold) = NaN;
    
    % I_ratioを計算（結果にはNaNが含まれる可能性がある）
    I_ratio = log(abs(Iup_processed ./ Idown_processed));
    I_ratio_list.(['ch' num2str(i)]) = I_ratio(1:datalen);
    
    % flow_velocityを計算（結果にはNaNが含まれる可能性がある）
    flow_velocity = I_ratio / K * cs *1e-3; 
    flow_list.(['ch' num2str(i)]) = flow_velocity(1:datalen);
end

%% --- プロットセクション ---

% % --- I_ratioのプロット ---
% f1 = figure('Name', ['I_ratio for Shot ' num2str(shotnum_orig)], 'NumberTitle', 'on');
% f1.WindowState = 'maximized';
% for i = 1:9
%     subplot(3, 3, i);
%     plot(I_ratio_list.time, I_ratio_list.(['ch' num2str(i)]));
%     grid on;
%     title(['Probe Head ', num2str(i)]);
%     xlabel('Time (us)');
%     ylabel('log(I_{up}/I_{down})');
%     xlim([470, 550]);
%     ylim([-5, 5]);
% end
% sgtitle(['I_{ratio} for Date: ' num2str(date) ', Shot: ' num2str(shotnum_orig)], 'FontSize', 14, 'FontWeight', 'bold');

% Iup, Idownのプロット
f2 = figure('Name', ['Processed Iup & Idown for Shot ' num2str(shotnum_orig)], 'NumberTitle', 'on');
f2.WindowState = 'maximized';
for i = 1:9
    subplot(3, 3, i);
    hold on;
    plot(Iup_raw_list.time, Iup_raw_list.(['ch' num2str(i)]), 'Color', 'magenta', 'DisplayName', 'I_{up} (Raw)');
    plot(Idown_raw_list.time, Idown_raw_list.(['ch' num2str(i)]), 'Color', 'cyan', 'DisplayName', 'I_{down} (Raw)');
    plot(Iup_list.time, Iup_list.(['ch' num2str(i)]), 'DisplayName', 'I_{up}', 'Color', 'red', 'LineWidth', 1);
    plot(Idown_list.time, Idown_list.(['ch' num2str(i)]), 'DisplayName', 'I_{down}', 'Color', 'blue', 'LineWidth', 1);
    
    yline(current_threshold_list(i), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Threshold');
    hold off;
    grid on;
    title(['Probe Head ', num2str(i)]);
    xlabel('Time (us)');
    ylabel('Current [V]');
    xlim([470, 550]);
    legend;
end
sgtitle(['Processed Currents for Date: ' num2str(date) ', Shot: ' num2str(shotnum_orig)], 'FontSize', 14, 'FontWeight', 'bold');

% Flow速度のプロット
f3 = figure('Name', ['Flow Velocity for Shot ' num2str(shotnum_orig)], 'NumberTitle', 'on');
f3.WindowState = 'maximized';
for i = 1:9
    subplot(3, 3, i);
    flow_to_plot = fillmissing(flow_list.(['ch' num2str(i)]), 'linear');
    plot(flow_list.time, flow_to_plot);
    grid on;
    title(['Probe Head ', num2str(i)]);
    xlabel('Time (us)');
    ylabel('Flow Velocity [km/s]');
    xlim([470, 550]);
    ylim([-100, 100]); % 必要に応じて調整してください
end
sgtitle(['Flow Velocity for Date: ' num2str(date) ', Shot: ' num2str(shotnum_orig)], 'FontSize', 14, 'FontWeight', 'bold');

disp('プロットが完了しました。');
