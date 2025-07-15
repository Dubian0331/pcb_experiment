% MP_Iratio関数を使って、指定範囲のshotnumでI_ratioファイルを一括生成するスクリプト

date = 240611;
start_num = 47;
end_num = 72;

for shotnum = start_num:end_num
    try
        I_ratio_list = MP_Iratio(date, shotnum);
        % f1 = figure('Name', ['I_ratio for Shot ' num2str(shotnum_orig)], 'NumberTitle', 'on');
        % f1.WindowState = 'maximized';
        % for i = 1:9
        %     subplot(3, 3, i);
        %     plot(I_ratio_list.time, I_ratio_list.(['ch' num2str(i)]));
        %     grid on;
        %     title(['Probe Head ', num2str(i)]);
        %     xlabel('Time (us)');
        %     ylabel('log(I_{up}/I_{down})');
        %     xlim([470, 550]);
        %     ylim([-5, 5]);
        % end
        % close(f1);
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
