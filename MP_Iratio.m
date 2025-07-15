%%
% ------------------------------------------------------------ %
% shotnumを指定して、上流と下流の電流の比をI_ratioというフォルダに格納
% ------------------------------------------------------------ %
function I_ratio_list = MP_Iratio(date, shotnum)
% clear all;
% close all;

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
fs = 1e7;
fc = 1e5;

% チャンネルごとにプロット
channel_list = [1, 3];
% channel_list = [2, 4];

I_ratio_list = table;
I_ratio_list.("time")=df.("time(us)")(1:datalen);
for i = 1:9
    ch_up_name = sprintf('ch%d', 4*(i-1) + channel_list(1));
    ch_down_name = sprintf('ch%d', 4*(i-1) + channel_list(2));

    % 生信号を代入
    Iup_raw = df.(ch_up_name);
    Idown_raw = df.(ch_down_name);

    % 負の値をNaNに置換
    Iup_raw(Iup_raw < 0) = NaN;
    Idown_raw(Idown_raw < 0) = NaN;

    % NaNを線形補間し、lowpassフィルタを適用
    Iup_filtered = lowpass(fillmissing(Iup_raw, 'linear'), fc, fs);
    Idown_filtered = lowpass(fillmissing(Idown_raw, 'linear'), fc, fs);

    % 5. フィルタ後のデータに、元のNaNの位置を復元
    Iup_processed = Iup_filtered;
    Idown_processed = Idown_filtered;
    Iup_processed(isnan(Iup_raw)) = NaN;
    Idown_processed(isnan(Idown_raw)) = NaN;

    % 6. 信号が小さすぎる区間を無効化するための閾値処理
    peak_current = max(max(Iup_processed, [], 'omitnan'), max(Idown_processed, [], 'omitnan'));
    if ~isempty(peak_current) && peak_current > 0
        current_threshold = peak_current * 0.01; % ピーク値の1%を閾値とする
        Iup_processed(Iup_processed < current_threshold) = NaN;
        Idown_processed(Idown_processed < current_threshold) = NaN;
    end
    I_ratio = log(abs(Iup_processed ./ Idown_processed));
    I_ratio = fillmissing(I_ratio, 'linear');
    I_ratio_list.(strcat("ch", num2str(i))) = I_ratio(1:datalen);
end

mkdir(strcat(filepath,sprintf('\\%d\\I_ratio', date)));
writetable(I_ratio_list, strcat(filepath,"\\",num2str(date),"\\I_ratio\\",shotnum,".csv"), 'WriteVariableNames', true);