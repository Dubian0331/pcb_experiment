%%
% % ------------------------------------------------------------ %
% % 240610 shotnum=35:47について、リレー回路の不具合により、CHを入れ替えたので修正
% % ------------------------------------------------------------ %
function MP_modify_Iratio(date, shotnum, src_dir, dst_dir)
    % shotnum: 数値またはゼロ埋め文字列
    % date: 数値または文字列
    % src_dir, dst_dir: フォルダパス
    if isnumeric(date)
        date = num2str(date);
    end
    if isnumeric(shotnum)
        if shotnum < 10
            shotstr = ['0', num2str(shotnum)];
        else
            shotstr = num2str(shotnum);
        end
    else
        shotstr = shotnum;
    end
    src_file = fullfile(src_dir, date, 'I_ratio', [shotstr, '.csv']);
    dst_folder = fullfile(dst_dir, date, 'I_ratio_modified');
    dst_file = fullfile(dst_folder, [shotstr, '.csv']);
    if ~isfile(src_file)
        warning('ファイルが見つかりません: %s', src_file);
        return;
    end
    data = readtable(src_file);
    G = data{:, 'ch6'};
    H = data{:, 'ch7'};
    I = data{:, 'ch8'};
    data{:, 'ch5'} = G;
    data{:, 'ch6'} = H;
    data{:, 'ch7'} = I;
    if ~exist(dst_folder, 'dir')
        mkdir(dst_folder);
    end
    writetable(data, dst_file);
    disp(['修正版を保存: ', dst_file]);
end


% % CSVファイルの読み込み
% filename = 'G:\My Drive\lab\lab_data\mach_probe_rawdata\240610\I_ratio\017.csv';
% data = readtable(filename);

% % G列、H列、I列を取得
% G = data{:, 'ch6'};
% H = data{:, 'ch7'};
% I = data{:, 'ch8'};

% % 列を置き換える
% data{:, "ch5"} = G;
% data{:, 'ch6'} = H;
% data{:, 'ch7'} = I;

% % 結果を新しいCSVファイルに保存
% new_filename = 'G:\My Drive\lab\lab_data\mach_probe_rawdata\240610\I_ratio_modified\017.csv';
% writetable(data, new_filename);

% %%
% function MP_Iratio(shotnum)
%     % データの読み込み
%     date = 240610;

%     if shotnum < 10
%         shotnum = ['00', num2str(shotnum)];
%     else 
%         shotnum = ['0', num2str(shotnum)];
%     end

%     % CSVデータの入っているフォルダの指定
%     filepath = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\MachProbe_data';
%     filename = strcat(filepath, '\', num2str(date), '\I_ratio\', shotnum, '.csv');

%     % データをCSVファイルから読み込む
%     data = readtable(filename);

%     % G列、H列、I列を取得
%     G = data{:, 'ch6'};
%     H = data{:, 'ch7'};
%     I = data{:, 'ch8'};

%     % 列を置き換える
%     data{:, 'ch5'} = G;
%     data{:, 'ch6'} = H;
%     data{:, 'ch7'} = I;

%     % 結果を新しいCSVファイルに保存
%     mkdir(strcat(filepath,'\', num2str(date), '\I_ratio_modified'));
%     new_filename = strcat(filepath, '\', num2str(date), '\I_ratio_modified\', shotnum, '.csv');
%     writetable(data, new_filename);
% end

% %%
% % 関数をループで呼び出す
% for i = 35:47
%     MP_Iratio(i)
% end
