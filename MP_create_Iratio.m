% MP_Iratio関数を使って、指定範囲のshotnumでI_ratioファイルを一括生成するスクリプト

start_num = 2;
end_num = 47;

for shotnum = start_num:end_num
    try
        MP_Iratio(shotnum);
    catch ME
        if strcmp(ME.identifier, 'MATLAB:readtable:OpenFailed') || contains(ME.message, 'No such file or directory')
            disp(['ファイルが見つかりません: shotnum = ', num2str(shotnum)]);
        else
            disp(['エラーが発生しました: shotnum = ', num2str(shotnum)]);
            disp(ME.message);
        end
        continue;
    end
end
