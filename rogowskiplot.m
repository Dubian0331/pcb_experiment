%%
% ------------------------------------------------------------ %
% ロゴスキー波形をプロット
% ------------------------------------------------------------ %
function [] = rogowskiplot(shotnum)
close all;

% データの読み込み
date = 240204;
axisfontsize = 12;

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

newTable = table(time, TF, PF1, PF2);
writetable(newTable, strcat(filepath,"\\",num2str(date),"\\rogowski\\figure\\",shotnum,".csv"));

f4 = figure('Name', strcat(shotnum,'_rogowski'), 'NumberTitle', 'on');

% ロゴスキーのデータを同じグラフにプロットしたくなった
% plot(time, TF, DisplayName='TF', LineWidth=1.5, LineStyle='-', color='black')
hold on;
% plot(time, PF1, DisplayName='PF1', LineStyle='-')
% plot(time, PF2, DisplayName='PF2', LineStyle='-')  
plot(time, FCPF1, DisplayName='FCPF1', LineStyle='-')   
plot(time, FCPF2, DisplayName='FCPF2', LineStyle='-')   
plot(time, FCTF1, DisplayName='FCTF1', LineStyle='-')   
plot(time, FCTF2, DisplayName='FCTF2', LineStyle='-')  
ylabel('Voltage [V]', 'FontSize', axisfontsize);

grid on;
legend('FontSize', 6, 'Location', 'northeast');
xrange = [0, 1000];
xlim(xrange);
ylim([-2, 2]);
set(gca, 'YColor', 'k'); % メモリの色を黒に設定

 
saveas(gcf, strcat(filepath,"\\",num2str(date),"\\rogowski\\figure\\",shotnum,".png"));