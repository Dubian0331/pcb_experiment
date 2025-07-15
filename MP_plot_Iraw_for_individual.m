%% ========================================================================
%  マッハプローブ生信号プロット用スクリプト
% =========================================================================

clearvars -except default_mach;
% close all;

% --- 1. 'default_mach'変数から設定を読み込み、ダイアログを表示 ---
if ~exist('default_mach','var'), default_mach=NaN; end
answer = customDialog(default_mach);
default_mach = answer;
date = str2double(answer.dateValue);
shotnum_val = str2double(answer.shotValue);
doSave = answer.saveFigureValue;

% --- 2. ファイルパスの構築とデータ処理 ---
if shotnum_val < 10
    shotnum_str_formatted = ['00',num2str(shotnum_val)];
else 
    shotnum_str_formatted = ['0',num2str(shotnum_val)];
end
% filename = sprintf('\\\\NIFS\\experiment\\results\\MachProbe\\%d\\ES_%d0%s.csv', date, date, shotnum_str_formatted);
filepath = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\MachProbe_data';
filename = strcat(filepath,'\',num2str(date),'\ES_',num2str(date),shotnum_str_formatted,'.csv');
fprintf('読み込みファイル: %s\n', filename);

try
    df = readtable(filename, 'Delimiter', ',', "VariableNamingRule", "preserve");
catch ME
    errordlg(sprintf('ファイルの読み込みに失敗しました。\n\nエラー: %s\nファイルパス: %s', ME.message, filename), 'ファイルエラー');
    return;
end

f1 = figure('Name', ['Shot #', num2str(shotnum_val), ' - All Channels'], 'NumberTitle', 'on');
f1.WindowState = 'maximized';
colors = ["red", "cyan", "blue", "magenta"];
delay = 4;
target_channel = [1, 3];

for i = 1:9
    subplot(3, 3, i);
    hold on;
    grid on;
    for j = target_channel
        channel_name = sprintf('ch%d', 4*(i-1) + j);
        label = sprintf('ch%d-%d', i, j);
        if ismember(channel_name, df.Properties.VariableNames)
            filtered_data = lowpass(df.(channel_name), 1e3, 10e6);
            plot(df.('time(us)') + delay, filtered_data, 'DisplayName', label, 'Color', colors(j), 'LineStyle', '-');
        end
    end
    xlabel('Time (us)');
    if any(i == [1, 4, 7])
        ylabel('Voltage (V)');
    end
    legend('FontSize', 7, 'Location', 'northeast');
    xlim([450, 600]);
end
sgtitle(['Date: ', num2str(date), ', Shot: ', num2str(shotnum_val)], 'FontSize', 14, 'FontWeight', 'bold');
disp('プロットが完了しました。');

% --- 保存 ---
if doSave
    saveDir = fullfile(filepath, num2str(date), 'figure');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    
    % 保存ファイル名を生成
    saveFilename = sprintf('%s_Iraw.png', shotnum_str_formatted);
    fullSavePath = fullfile(saveDir, saveFilename);
    
    % グラフを保存
    saveas(f1, fullSavePath);
    fprintf('グラフを保存しました: %s\n', fullSavePath);
end

function answer = customDialog(default)
    if ~isfield(default,'dateValue')
        defaultDate = ''; 
    else
        defaultDate = default.dateValue;
    end
    if ~isfield(default,'shotValue')
        defaultShot = ''; 
    else
        defaultShot = default.shotValue;
    end
    if ~isfield(default,'saveFigureValue')
        defaultSave = false; 
    else
        defaultSave = default.saveFigureValue;
    end

    % ダイアログボックスを作成
    d = dialog('Position', [300, 300, 400, 240], 'Name', 'Custom Dialog');

    % 'date'テキスト入力フィールドのラベル
    uicontrol('Parent', d, ...
              'Style', 'text', ...
              'Position', [20, 180, 100, 20], ...
              'String', 'Date:');
    
    % 'date'テキスト入力フィールド
    dateField = uicontrol('Parent', d, ...
                          'Style', 'edit', ...
                          'Position', [130, 180, 200, 25], ...
                          'String', defaultDate);

    % 'shot'テキスト入力フィールドのラベル
    uicontrol('Parent', d, ...
              'Style', 'text', ...
              'Position', [20, 150, 100, 20], ...
              'String', 'Shot:');
    
    % 'shot'テキスト入力フィールド
    shotField = uicontrol('Parent', d, ...
                          'Style', 'edit', ...
                          'Position', [130, 150, 200, 25], ...
                          'String', defaultShot);

    % 'save'チェックボックス
    saveCheckbox = uicontrol('Parent', d, 'Style', 'checkbox', 'Position', [130, 100, 200, 25], 'String', 'Save', 'Value', defaultSave);

    % OKボタン
    btn = uicontrol('Parent', d, ...
                    'Position', [150, 40, 70, 25], ...
                    'String', 'OK', ...
                    'Callback', @btnCallback);

    % ボタンのコールバック関数
    function btnCallback(~, ~)
        % 各フィールドの値を取得
        answer.dateValue = get(dateField, 'String');
        answer.shotValue = get(shotField, 'String');
        answer.saveFigureValue = get(saveCheckbox, 'Value');
        % ダイアログボックスを閉じる
        delete(d);
    end

    % ダイアログボックスが閉じられるまで待機
    uiwait(d);

    % 関数の返り値として出力
    if isvalid(d)
        delete(d);
    end
end