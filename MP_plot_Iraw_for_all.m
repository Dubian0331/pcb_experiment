%% ========================================================================
%  マッハプローブ生信号プロット用スクリプト
% =========================================================================
clearvars -except default_mach;
close all;

% --- 1. 'default_mach'変数から設定を読み込み、ダイアログを表示 ---
if ~exist('default_mach', 'var') || ~isstruct(default_mach)
    default_mach = struct(); % 存在しない場合は空の構造体として作成
end

answer = customDialog_batch(default_mach);

% ダイアログがキャンセルされた場合はスクリプトを終了
if isempty(answer)
    disp('処理がキャンセルされました。');
    return;
end

% 今回の入力値を次回のために記憶
default_mach = answer;

% ダイアログからの値を取得
date = str2double(answer.dateValue);
start_shot = str2double(answer.startShotValue);
end_shot = str2double(answer.endShotValue);


%% --- 2. 指定範囲のショットをループ処理 ---
fprintf('バッチ処理を開始します: Date %d, Shot %d から %d まで\n', date, start_shot, end_shot);

for shotnum_val = start_shot:end_shot
    
    fprintf('--> 処理中: Shot #%d...', shotnum_val);
    
    % --- ファイルパスの構築 ---
    if shotnum_val < 10
        shotnum_str_formatted = ['00',num2str(shotnum_val)];
    elseif shotnum_val < 100
        shotnum_str_formatted = ['0',num2str(shotnum_val)];
    else
        shotnum_str_formatted = num2str(shotnum_val);
    end
    
    filepath = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\MachProbe_data';
    filename = strcat(filepath,'\',num2str(date),'\ES_',num2str(date),shotnum_str_formatted,'.csv');

    % --- データの読み込みとプロット ---
    try
        df = readtable(filename, 'Delimiter', ',', "VariableNamingRule", "preserve");
    catch
        fprintf(' スキップ (ファイルが見つかりません)\n');
        continue; % 次のループへ
    end
    
    % Figureを画面に表示しない
    f1 = figure('Name', ['Shot #', num2str(shotnum_val)], 'NumberTitle', 'off', 'Visible', 'off');
    
    colors = ["red", "cyan", "blue", "magenta"];
    delay = 0;
    target_channel = [1, 3];

    for i = 1:9
        subplot(3, 3, i);
        hold on;
        grid on;
        for j = target_channel
            channel_name = sprintf('ch%d', 4*(i-1) + j);
            label = sprintf('ch%d-%d', i, j);
            if ismember(channel_name, df.Properties.VariableNames)
                raw_signal = df.(channel_name);
                clean_signal = fillmissing(raw_signal, 'linear');
                filtered_data = lowpass(clean_signal, 1e6, 10e6);
                % filtered_data = lowpass(df.(channel_name), 1e6, 10e6);
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
    
    % --- グラフの保存 ---
    saveDir = fullfile(filepath, num2str(date), 'figure');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    saveFilename = sprintf('%s_Iraw.png', shotnum_str_formatted);
    fullSavePath = fullfile(saveDir, saveFilename);
    saveas(f1, fullSavePath);
    
    close(f1);
    
    fprintf(' 保存完了\n');
    
end

disp('全ての処理が完了しました。');


% --- ダイアログ関数 ---
function answer = customDialog_batch(default)
    answer = [];

    % --- デフォルト値の準備 ---
    if ~isstruct(default) || ~isfield(default, 'dateValue')
        defaultDate = '240611'; 
    else
        defaultDate = default.dateValue;
    end
    
    % 前回の終了ショットの次を開始ショットのデフォルトにする
    if ~isstruct(default) || ~isfield(default, 'endShotValue') || isempty(default.endShotValue)
        defaultStartShot = ''; 
    else
        prev_end_shot = str2double(default.endShotValue);
        defaultStartShot = num2str(prev_end_shot + 1);
    end
    
    defaultEndShot = ''; % 終了ショットは毎回入力

    % ダイアログボックスを作成
    d = dialog('Position', [300, 300, 400, 240], 'Name', 'バッチ処理設定');
    
    uicontrol('Parent', d, 'Style', 'text', 'Position', [20, 180, 100, 20], 'String', 'Date:');
    dateField = uicontrol('Parent', d, 'Style', 'edit', 'Position', [130, 180, 200, 25], 'String', defaultDate);
    
    uicontrol('Parent', d, 'Style', 'text', 'Position', [20, 140, 100, 20], 'String', '開始ショット:');
    startShotField = uicontrol('Parent', d, 'Style', 'edit', 'Position', [130, 140, 200, 25], 'String', defaultStartShot);
    
    uicontrol('Parent', d, 'Style', 'text', 'Position', [20, 100, 100, 20], 'String', '終了ショット:');
    endShotField = uicontrol('Parent', d, 'Style', 'edit', 'Position', [130, 100, 200, 25], 'String', defaultEndShot);
    
    uicontrol('Parent', d, 'Position', [150, 40, 70, 25], 'String', '実行', 'Callback', @btnCallback);

    function btnCallback(~, ~)
        answer.dateValue = get(dateField, 'String');
        answer.startShotValue = get(startShotField, 'String');
        answer.endShotValue = get(endShotField, 'String');
        
        % 入力が空でないかチェック
        if isempty(answer.dateValue) || isempty(answer.startShotValue) || isempty(answer.endShotValue)
            warndlg('全てのフィールドを入力してください。', '入力エラー');
            answer = []; % answerを空に戻す
        end
        delete(d);
    end
    
    uiwait(d);
    if ishandle(d)
        delete(d);
    end
end