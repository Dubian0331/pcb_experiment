clear all;
close all;

titlefontsize = 20;
axisfontsize = 12;
delay = 4; % トリガータイミングがずれてそう？


% ------------------------------------------------------------ %
date = 240204; 
r = 15; 
firstshot = 6;

% % csvデータの入っているフォルダの指定
filepath = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data';

% Excelファイルを読み込む
filename = strcat(filepath, "\\", num2str(date), '\\flow_list_r=', num2str(r), '.xlsx');

% シート名のリスト
sheet_names = {'06', '07', '08', '09', '10'};

% 列のリスト
column_names = {'ch1', 'ch2', 'ch3', 'ch4', 'ch5', 'ch6', 'ch7', 'ch8', 'ch9'};

% 新しいExcelファイル名
output_filename = strcat(filepath, '\\', num2str(date), '\\data_r=', num2str(r), '.xlsx');

% 新しいExcelファイルにデータを書き込む
data1 = 
data2 = 
for c = 1:length(column_names)
    data_table = table;
    for s = 1:length(sheet_names)
        % シート名を指定してデータを読み込む
        data = readtable(filename, 'Sheet', sheet_names{s}, );

        if s == 1
            % time列を各シートの一番左に追加
            data = [data(:, {'time'}), data(:, column_names{c})];
            % 変数名を変更
            data.Properties.VariableNames = {'time', num2str(s)};
        else
            data = data(:, column_names{c});
            data.Properties.VariableNames = {num2str(s)};
        end
        % 
        % % テーブルをセルに変換
        % data_cell{s} = data;

        

        % テーブルにデータを追加
        if isempty(data_table)
            data_table = data;
        else
            data_table = [data_table, data];
        end
    end

    

    % 変数名を追加
    variable_names = [{'time'}, {'1'}, {'2'}, {'3'}, {'4'}, {'5'}];
    % data_table.Properties.VariableNames = variable_names;

    % テーブルをExcelファイルに書き込む
    writetable(data_table, output_filename, 'Sheet', column_names{c}, 'WriteVariableNames', true);
end

%%
for i = 1:9
    filename = strcat(filepath, '\\', num2str(date), '\\data_r=', num2str(r), '.xlsx')
    df = readtable(filename, 'Sheet', strcat("ch", num2str(i)));
    for j = 1:5
        a += df.(num2str(j));
end


%% 
