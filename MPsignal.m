%%

function [] = MP_signal(shotnum)

close all;

axisfontsize = 12;
delay = 4; % トリガータイミングがずれてそう？
% ------------------------------------------------------------ %
% データの読み込み
date = 240204;

if shotnum < 10
    shotnum = sprintf('%02d', shotnum);
else 
    shotnum = num2str(shotnum);
end

% csvデータの入っているフォルダの指定
filepath = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data';
filename = strcat(filepath,sprintf('\\%d\\ES_%d0%s.csv', date, date, shotnum));
% filename = sprintf('\\\\NIFS\\experiment\\results\\MachProbe\\%d\\ES_%d0%s.csv', date, date, shotnum);

df = readtable(filename, 'Delimiter', ',',"VariableNamingRule","preserve");

% 生信号をプロット

% プロットの準備
f1 = figure('Name', strcat(shotnum,'_All Channels'), 'NumberTitle', 'on');
colors = ["red", "cyan", "blue", "magenta"];

% チャンネルごとにプロット
f1.WindowState = 'maximized';
channel_list = [];
channel_list = [channel_list, 1];
% channel_list = [channel_list, 2];
channel_list = [channel_list, 3];
% channel_list = [channel_list, 4];

for i = 1:9
    subplot(3, 3, i);
    hold on;
    for j = channel_list
        channel_name = sprintf('ch%d', 4*(i-1) + j);
        label = sprintf('ch%d-%d', i, j);
        % filtered_data = lowpass(df.(channel_name), 1e6, 10e6);
        plot(df.('time(us)')+delay, lowpass(df.(channel_name), 1e6, 10e6), 'DisplayName', label, 'Color', colors(j), LineStyle='-');
    end

    % ラベルや凡例の設定
    xlabel('Time (us)', 'FontSize', axisfontsize);
    if any(i == [1, 4, 7])
        ylabel('Voltage (V)', 'FontSize', axisfontsize);
    end
    grid on;
    legend('FontSize', 6, 'Location', 'northeast');
    % xrange = [450, 600];
    xrange = [475, 525];
    % xrange = [390, 405];
    % xrange = [0, 1000];
    xlim(xrange);
    yrange_left = [-0.5, 0.5];
    ylim(yrange_left);
    set(gca, 'YColor', 'k'); % メモリの色を黒に設定
end

mkdir(strcat(filepath,sprintf('\\%d\\figure', date)));
saveas(gcf, strcat(filepath, sprintf('\\%d\\figure\\%s_all_%d-%dus.png', date, shotnum, xrange(1), xrange(2))), 'png');