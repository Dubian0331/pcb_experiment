titlefontsize = 20;
axisfontsize = 12;

N = 200;
R = 200; 
delay = 0; 

% データの読み込み
date = 241229;
shotnum = 31;
if shotnum < 10
    shotnum = sprintf('%02d', shotnum);
else 
    shotnum = num2str(shotnum);
end

% filename = sprintf('//NIFS/experiment/results/MachProbe/%d/ES_%d0%s.csv', date, date, shotnum);
% filepath = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data';
% filename = strcat(filepath,sprintf('\\%d\\ES_%d0%s.csv', date, date, shotnum));
filename = sprintf('\\\\NIFS\\experiment\\results\\PotentialProbe\\%d\\PP_%d0%s.csv', date, date, shotnum);
df = readtable(filename, 'Delimiter', ',',"VariableNamingRule","preserve");

%% 

% prompt = 'Enter the channel number (1-9): ';
% dlgtitle = 'Channel Selection';
% dims = [1 50];
% definput = {''};
% answer = inputdlg(prompt, dlgtitle, dims, definput);
% 生信号をプロット

% プロットの準備
f1 = figure('Name', strcat(shotnum,'_PP'), 'NumberTitle', 'on');
colors = ["red", "cyan", "blue", "magenta"];

% チャンネルごとにプロット
f1.WindowState = 'maximized';
channel_list = [];
channel_list = [channel_list, 1];
channel_list = [channel_list, 2];
channel_list = [channel_list, 3];
channel_list = [channel_list, 4];

for i = 1:8
    subplot(3, 3, i);
    hold on;
    channel_name = sprintf('ch%d', i);
    label = channel_name;
    % yyaxis left;
    % plot(df.('time(us)')+delay, smoothdata(df.(channel_name)), 'DisplayName', label, 'Color', colors(j), LineStyle='-');
    filtered_data = lowpass(df.(channel_name), 1e6, 1e7);
    filtered_data2 = smoothdata(filtered_data);
    plot(df.('time(us)')+delay, filtered_data2, 'DisplayName', label, LineStyle='-');

    % ラベルや凡例の設定
    xlabel('Time (us)', 'FontSize', axisfontsize);
    if any(i == [1, 4, 7])
        ylabel('Voltage (V)', 'FontSize', axisfontsize);
    end
    grid on;
    legend('FontSize', 6, 'Location', 'northeast');
    % xrange = [300, 600];
    xrange = [400, 525];
    % xrange = [390, 405];
    % xrange = [0, 1000];
    xlim(xrange);
    % yrange_left = [-0.5, 0.5];
    yrange_left = [-3, 3];
    ylim(yrange_left);
    set(gca, 'YColor', 'k'); % メモリの色を黒に設定
end