%%
clear all;
% close all;

titlefontsize = 12;
axisfontsize = 10;
N = 200;
R = 200;
delay = 4; % トリガータイミングがずれてそう？
% ------------------------------------------------------------ %
% データの読み込み
date = 240611;


prompt = 'Enter the shot number: ';
dlgtitle = 'Shot Number';
dims = [1 30];
definput = {''};
shotnum = inputdlg(prompt, dlgtitle, dims, definput);
shotnum = str2double(shotnum{1});

if shotnum < 10
    shotnum = sprintf('%02d', shotnum);
else 
    shotnum = num2str(shotnum);
end

% csvデータの入っているフォルダの指定
% filepath = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data';
% filename = strcat(filepath,sprintf('\\%d\\ES_%d0%s.csv', date, date, shotnum));
filename = sprintf('\\\\NIFS\\experiment\\results\\MachProbe\\%d\\ES_%d0%s.csv', date, date, shotnum);
% filename = sprintf('\\\\NIFS\\experiment\\results\\ElectroStaticProbe\\%d\\ES_%d0%s.csv', date, date, shotnum);
% filename = sprintf('\\\\NIFS\\experiment\\results\\PotentialProbe\\%d\\PP_%d0%s.csv', date, date, shotnum);

df = readtable(filename, 'Delimiter', ',',"VariableNamingRule","preserve");
% 生信号をプロット

% プロットの準備
f1 = figure('Name', strcat(shotnum,'_All Channels'), 'NumberTitle', 'on');
colors = ["red", "cyan", "blue", "magenta"];

% チャンネルごとにプロット
f1.WindowState = 'maximized';
channel_list = [];
channel_list = [channel_list, 1];
channel_list = [channel_list, 2];
channel_list = [channel_list, 3];
channel_list = [channel_list, 4];

for i = 1:9
    subplot(3, 3, i);
    hold on;
    for j = channel_list
        channel_name = sprintf('ch%d', 4*(i-1) + j);
        label = sprintf('ch%d-%d', i, j);
        % yyaxis left;
        % plot(df.('time(us)')+delay, smoothdata(df.(channel_name)), 'DisplayName', label, 'Color', colors(j), LineStyle='-');
        filtered_data = lowpass(df.(channel_name), 1e6, 10e6);
        plot(df.('time(us)')+delay, lowpass(df.(channel_name), 1e6, 10e6), 'DisplayName', label, 'Color', colors(j), LineStyle='-');
    end

    % ラベルや凡例の設定
    xlabel('Time (us)', 'FontSize', axisfontsize);
    if any(i == [1, 4, 7])
        ylabel('Voltage (V)', 'FontSize', axisfontsize);
    end
    grid on;
    legend('FontSize', 6, 'Location', 'northeast');
    xrange = [450, 600];
    % xrange = [475, 525];
    % xrange = [390, 405];
    % xrange = [0, 1000];
    xlim(xrange);
    yrange_left = [-0.5, 0.5];
    % ylim(yrange_left);
    set(gca, 'YColor', 'k'); % メモリの色を黒に設定
end

% グラフの保存
% mkdir(strcat(filepath,sprintf('\\%d\\figure', date)));
% saveas(gcf, strcat(filepath, sprintf('\\%d\\figure\\%s_all_%d-%dus.png', date, shotnum, xrange(1), xrange(2))), 'png');

%% 
% フローz-側（シールドルーム側）から見て時計回りが+のフロー

%　フローをプロット
K = 2.5; % 比例定数
k = 1.38e-23;    % [J/K]
q = 1.60e-19;  % [C]
Te = 5 * q / k; % [K]　% T_e：電子温度→青山さんのp.35を参照 5~8[eV]
mi = 1.672e-27;   % [kg]
cs = sqrt(2 * k * Te / mi);

% プロットの準備
f3 = figure('Name', strcat(shotnum,'_flow'), 'NumberTitle', 'on');
colors = ["cyan", "red", "magenta", "blue"];
datalen = 10001;

% チャンネルごとにプロット
f3.WindowState = 'maximized';
channel_list = [1, 3];
% channel_list = [2, 4];

Iup = [];
Idown = [];
Zup = channel_list(1);
Zdown = channel_list(2);
flow_list = table;
I_ratio_list = table;
flow_list.("time")=df.("time(us)")(1:datalen);
I_ratio_list.("time")=df.("time(us)")(1:datalen);
for i = [1:5, 7:8] 
    subplot(3, 3, i);
    hold on;
    for j = channel_list
        channel_name = sprintf('ch%d', 4*(i-1) + j);
        if j == channel_list(1)
            % Iup = df.(channel_name);
            Iup = lowpass(df.(channel_name), 1e6, 10e6);
        else
            % Idown = df.(channel_name);
            Idown = lowpass(df.(channel_name), 1e6, 10e6);
        end
    end


    Flow = log(Iup ./ Idown) / K * cs * 1e-3;  % [km/s]
    valid_indices = imag(Flow) == 0;
    valid_flow = Flow(valid_indices);
    valid_time = df.('time(us)')(valid_indices);
    % plot(df.('time(us)'), Flow)
    plot(valid_time+delay, valid_flow)
    % hold on;
    for l = 1:numel(Flow)
        if imag(Flow(l)) ~= 0
            Flow(l) = 0;
        end
    end
    flow_list.(strcat("ch", num2str(i))) = Flow(1:datalen);
    I_ratio = log(Iup ./ Idown);
    for l = 1:numel(I_ratio)
        if imag(I_ratio(l)) ~= 0
            I_ratio(l) = 0;
        end
    end
    I_ratio_list.(strcat("ch", num2str(i))) = I_ratio(1:datalen);
    % ラベルや凡例の設定
    xlabel('Time (us)', 'FontSize', 12);
    if any(i == [1, 4, 7])
        ylabel('Vz (km/s)', 'FontSize', 12);
    end
    grid on;
    xrange = [484, 499];
    xlim(xrange);
    yrange = [-40, 40];
    ylim(yrange);
    title(strcat('ch', num2str(i)), 'Fontsize',titlefontsize);
    % set(f3, 'Name', strcat(shotnum,'_flow_z-=',num2str(Zup),'_z+=',num2str(Zdown)));
    set(f3, 'Name', strcat(shotnum,'_flow_-=',num2str(Zup),'_+=',num2str(Zdown)));
end
 
writetable(I_ratio_list, strcat(filepath,"\\",num2str(date),"\\I_ratio\\",shotnum,".csv"), 'WriteVariableNames', true);

% グラフの保存
% mkdir(strcat(filepath,sprintf('\\%d', date)));
% saveas(gcf, strcat(filepath, sprintf('\\%d/figure\\%s_flow_%d-%dus.png', date, shotnum, xrange(1), xrange(2))), 'png');

%%
for i = [6:35 37:41]
    MP_Iratio(i)
end

%%
for i = [6:35 37:41]
    rogowskiplot(i)
end

%%
for i = [6:35 37:41]
    MPsignal(i)
end

%%
% ロゴスキーのグラフをプロット
% データの読み込み
date = 240426;
titlefontsize = 20;
axisfontsize = 12;

prompt = 'Enter the shot number: ';
dlgtitle = 'Shot Number';
dims = [1 30];
definput = {''};
shotnum = inputdlg(prompt, dlgtitle, dims, definput);
shotnum = str2double(shotnum{1});

if shotnum < 10
    shotnum = sprintf('%02d', shotnum);
else 
    shotnum = num2str(shotnum);
end

% csvデータの入っているフォルダの指定
filepath = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data';

rogowskipath = strcat('\\rogowski\\', num2str(date), '0', shotnum, '.rgw');
rogowski = strcat(filepath,sprintf('\\%d', date), rogowskipath);
df2 = readtable(rogowski, 'Delimiter', '\t',"VariableNamingRule","preserve", "FileType","text");
time = df2.("time");
TF = df2.("ch1");
PF1 = df2.("ch3");
PF2 = df2.("ch4");
FCTF1 = df2.("ch5");
FCPF1 = df2.("ch7");
FCPF2 = df2.("ch8");
FCTF2 = df2.("ch9");

f4 = figure('Name', strcat(shotnum,'_rogowski'), 'NumberTitle', 'on');
plot(time, TF, DisplayName='TF', LineWidth=1.5, LineStyle='-', color='black')
hold on;
% plot(time, PF1, DisplayName='PF1', LineStyle='-')
% plot(time, PF2, DisplayName='PF2', LineStyle='-')  
% plot(time, FCPF1, DisplayName='FCPF1', LineStyle='-')   
% plot(time, FCPF2, DisplayName='FCPF2', LineStyle='-')   
% plot(time, FCTF1, DisplayName='FCTF1', LineStyle='-')   
% plot(time, FCTF2, DisplayName='FCTF2', LineStyle='-')  
ylabel('Voltage [V]', 'FontSize', axisfontsize);
grid on;
legend('FontSize', 6, 'Location', 'northeast');
xrange = [0, 1000];
xlim(xrange);
ylim([-2, 2]);
set(gca, 'YColor', 'k'); % メモリの色を黒に設定


%%
% プロットの準備
f1 = figure('Name', strcat(shotnum,'_All Channels'), 'NumberTitle', 'on');
colors = ["red", "cyan", "blue", "magenta"];

% チャンネルごとにプロット
f1.WindowState = 'maximized';

for i = 1:22
    channel_name = sprintf('ch%d', i);
    % yyaxis left;
    % plot(df.('time(us)')+delay, smoothdata(df.(channel_name)), 'DisplayName', label, 'Color', colors(j), LineStyle='-');
    % filtered_data = lowpass(df.(channel_name), 1e6, 10e6);
    plot(df.('time(us)')+delay, df.(channel_name), LineStyle='-', DisplayName=channel_name);
    % ラベルや凡例の設定
    xlabel('Time (us)', 'FontSize', axisfontsize);
    hold on;
    grid on;
    legend('FontSize', 6, 'Location', 'northeast');
    xrange = [450, 600];
    % xrange = [475, 525];
    % xrange = [390, 405];
    xrange = [0, 1000];
    xlim(xrange);
    yrange_left = [-0.5, 0.5];
    ylim(yrange_left);
    % set(gca, 'YColor', 'k'); % メモリの色を黒に設定
end