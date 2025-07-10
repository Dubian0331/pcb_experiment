%%
% ------------------------------------------------------------ %
% shotnumを指定して、上流と下流の電流の比をI_ratioというフォルダに格納
% ------------------------------------------------------------ %
function MP_Iratio(shotnum)
% clear all;
% close all;

% データの読み込み
date = 240610;

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

% プロットの準備
datalen = 10001;

% チャンネルごとにプロット
channel_list = [1, 3];
% channel_list = [2, 4];

Iup = [];
Idown = [];
I_ratio_list = table;
I_ratio_list.("time")=df.("time(us)")(1:datalen);
for i = 1:9
    for j = channel_list
        channel_name = sprintf('ch%d', 4*(i-1) + j);
        if j == channel_list(1)
            Iup = lowpass(df.(channel_name), 1e6, 10e6);
        else
            Idown = lowpass(df.(channel_name), 1e6, 10e6);
        end
    end

    I_ratio = log(Iup ./ Idown);
    for l = 1:numel(I_ratio)
        if imag(I_ratio(l)) ~= 0
            I_ratio(l) = 0;
        end
    end
    I_ratio_list.(strcat("ch", num2str(i))) = I_ratio(1:datalen);
end

mkdir(strcat(filepath,sprintf('\\%d\\I_ratio', date)));
writetable(I_ratio_list, strcat(filepath,"\\",num2str(date),"\\I_ratio\\",shotnum,".csv"), 'WriteVariableNames', true);