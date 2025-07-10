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

filename = sprintf('//NIFS/experiment/results/ElectroStaticProbe/%d/ES_%d0%s.csv', date, date, shotnum);
% filepath = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data';
filename = strcat(filepath,sprintf('\\%d\\ES_%d0%s.csv', date, date, shotnum));
df = readtable(filename, 'Delimiter', ',',"VariableNamingRule","preserve");

%% 

prompt = 'Enter the channel number (1-21): ';
dlgtitle = 'Channel Selection';
dims = [1 50];
definput = {''};
answer = inputdlg(prompt, dlgtitle, dims, definput);

% 入力されたチャンネル番号を取得
CHANNEL = str2double(answer{1});

% 入力値が1から9の範囲であるかを確認する
if isnan(CHANNEL) || CHANNEL < 1 || CHANNEL > 21 || rem(CHANNEL, 1) ~= 0
    error('Invalid channel number. Please enter a number between 1 and 9.');
end

% プロットの準備
f2 = figure('Name', strcat('CH', num2str(CHANNEL)), 'NumberTitle', 'on');
colors = ["cyan", "red", "magenta", "blue"];

% チャンネルごとにプロット
f2.WindowState = 'maximized';
ChannelName = sprintf('ch%d', CHANNEL);
Label = sprintf('ch%d', CHANNEL);
plot(df.('time(us)'), df.(ChannelName), 'DisplayName', Label);

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

