titlefontsize = 20;
axisfontsize = 12;

N = 200;
R = 200; 

% データの読み込み
date = 240128;
shotnum = 1;
if shotnum < 10
    shotnum = sprintf('%02d', shotnum);
else 
    shotnum = num2str(shotnum);
end

% filename = sprintf('//NIFS/experiment/results/MachProbe/%d/ES_%d0%s.csv', date, date, shotnum);
filepath = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data';
filename = strcat(filepath,sprintf('\\%d\\ES_%d0%s.csv', date, date, shotnum));
df = readtable(filename, 'Delimiter', ',',"VariableNamingRule","preserve");

%% 

prompt = 'Enter the channel number (1-9): ';
dlgtitle = 'Channel Selection';
dims = [1 50];
definput = {''};
answer = inputdlg(prompt, dlgtitle, dims, definput);

% 入力されたチャンネル番号を取得
CHANNEL = str2double(answer{1});

% 入力値が1から9の範囲であるかを確認する
if isnan(CHANNEL) || CHANNEL < 1 || CHANNEL > 9 || rem(CHANNEL, 1) ~= 0
    error('Invalid channel number. Please enter a number between 1 and 9.');
end

% プロットの準備
f2 = figure('Name', strcat('CH', num2str(CHANNEL)), 'NumberTitle', 'on');
colors = ["cyan", "red", "magenta", "blue"];
channel_list = [];
channel_list = [channel_list, 1];
channel_list = [channel_list, 2];
channel_list = [channel_list, 3];
channel_list = [channel_list, 4];

% チャンネルごとにプロット
f2.WindowState = 'maximized';
for i = channel_list
    ChannelName = sprintf('ch%d', 4*(CHANNEL-1)+i);
    Label = sprintf('ch%d-%d', CHANNEL, i);
    plot(df.('time(us)'), df.(ChannelName), 'DisplayName', Label, 'Color', colors{i});
    hold on;
end

xlabel('time [us]', 'FontSize', axisfontsize);
ylabel('Voltage [V]', 'FontSize', axisfontsize);

grid on;
legend('Location', 'northeast', 'FontSize', 6);
xrange = [430, 550];
% xrange = [0, 2000];
xlim(xrange);
yrange = [-2, 2];
ylim(yrange);
hold off;

% グラフの保存
% "データの入っているフォルダ/日付/figure"に保存される
% mkdir(strcat(filepath, sprintf('\\%d\\figure', date)));
% saveas(gcf, strcat(filepath, sprintf('\\%d\\figure\\%s_channel%d_%d-%dus.png', date, shotnum, CHANNEL, xrange(1), xrange(2))), 'png');  % グラフの保存

