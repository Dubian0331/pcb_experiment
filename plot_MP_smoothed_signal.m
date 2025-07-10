%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% スムージングと移動平均あり
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all

titlefontsize = 30;
axisfontsize = 20;

N = 200;
R = 200;
delay_time = 0;

% データの読み込み
date = 241230;
shotnum = 6

save_filepath = "G:\My Drive\lab\lab_data\mach_probe_rawdata";

filepath = 'G:\My Drive\lab\lab_data\mach_probe_rawdata';
% filepath = '\\NIFS\experiment\results\MachProbe';


if shotnum < 10
    shotnum_str = ['0', num2str(shotnum)];  % 1桁の場合はゼロパディング
    filename = strcat(filepath, '\', num2str(date), '\ES_', num2str(date), '0', shotnum_str, '.csv');
elseif shotnum < 100
    shotnum_str = num2str(shotnum, '%02d');  % 2桁に整形（代替法）
    filename = strcat(filepath, '\', num2str(date), '\ES_', num2str(date), '0', shotnum_str, '.csv');
else
    shotnum_str = num2str(shotnum);  % そのまま文字列化
    filename = strcat(filepath, '\', num2str(date), '\ES_', num2str(date), shotnum_str, '.csv');
end


df = readtable(filename, 'Delimiter', ',', "VariableNamingRule", "preserve");

%%%%%%%%%%%%%%%%% scale dataの読み込み %%%%%%%%%%%%%%%%%%
% スケール因子データの読み込み
scale_filepath = "G:\My Drive\lab\lab_data\mach_probe\analysis\scale_factors.xlsx"; % スケール因子が記載されたファイル
% 最新の年月日のシート名を取得
[~, sheet_names] = xlsfinfo(scale_filepath); % Excelファイル内の全シート名を取得
date_numbers = cellfun(@(x) str2double(x), sheet_names, 'UniformOutput', false); % 数値形式に変換可能なシート名を探す
date_numbers = cell2mat(date_numbers(~cellfun(@isempty, date_numbers))); % 有効な数値だけを抽出
[~, latest_index] = max(date_numbers); % 最新の日付を持つシートのインデックスを取得
latest_sheet = sheet_names{latest_index}; % 最新のシート名
% 最新のシートからデータを読み込む
scale_table = readtable(scale_filepath, 'Sheet', latest_sheet); % 最新のシートからスケール因子データを読み込み

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 信号処理パラメータ
Fs = 1e8; % サンプリング周波数 (1 MHzと仮定)
Fc = 1e6; % カットオフ周波数 (10 kHzと仮定)
windowSize = 100; % スムージング窓のサイズ (移動平均用)

% プロットの準備
f1 = figure('Name', strcat(shotnum,'_All Channels Processed'), 'NumberTitle', 'on');
colors = ["red", "cyan", "blue", "magenta"];

% チャンネルごとにプロット
f1.WindowState = 'maximized';
channel_list = [];
channel_list = [channel_list, 1];
channel_list = [channel_list, 3];

for i = 1:9
    subplot(3, 3, i);
    hold on;
    for j = channel_list
        channel_name = sprintf('ch%d', 4*(i-1) + j);
        label = sprintf('ch%d-%d', i, j);

        % スケール因子を取得
        scale_factor = scale_table.Scale(strcmp(scale_table.Channel, channel_name)); % 該当するスケール因子を取得

        % 元信号の取得
        time = df.('time(us)') + delay_time;
        signal = df.(channel_name) * scale_factor;

        % ローパスフィルタ処理
        filtered_signal = lowpass(signal, Fc, Fs);

        % スムージング処理 (移動平均フィルタ)
        smoothed_signal = movmean(filtered_signal, windowSize);

        % プロット
        plot(time, smoothed_signal, 'DisplayName', label, 'Color', colors(j));
    end

    % ラベルや凡例の設定
    xlabel('Time (us)', 'FontSize', axisfontsize);
    if any(i == [1, 4, 7])
        ylabel('Voltage (V)', 'FontSize', axisfontsize);
    end
    grid on;
    legend('FontSize', 15, 'Location', 'northeast');
    xrange = [450, 500];
    xlim(xrange);
    yrange = [-0.01, 0.05];
    ylim(yrange);
end

% フィギュア全体のタイトルを設定
sgtitle(sprintf('%d Shot %d', date, shotnum), 'FontSize', 20);

% グラフの保存
if not(exist(strcat(filepath,sprintf('/%d/figure', date))))
    mkdir(strcat(filepath,sprintf('/%d/figure', date)));
end
saveas(gcf, strcat(filepath, sprintf('/%d/figure/%s_all_%d-%dus_processed.png', date, shotnum_str, xrange(1), xrange(2))), 'png');
disp("saved processed figure")
