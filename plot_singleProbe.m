% singleProbeフォルダ内の全CSVファイルをプロットするスクリプト

function plot_singleProbe_main
    persistent defaultInput
    if isempty(defaultInput)
        defaultInput.dateValue = '';
        defaultInput.shotValue = '';
    end
    answer = customDialog(defaultInput);
    if isempty(answer)
        return;
    end
    defaultInput = answer;
    folder_date = answer.dateValue;
    file_num = str2double(answer.shotValue);
    if isnan(file_num)
        errordlg('Shot番号は整数で入力してください。', '入力エラー');
        return;
    end
    file_name = sprintf('%03d.csv', file_num);
    folder_path = fullfile('singleProbe', folder_date);
    file_path = fullfile(folder_path, file_name);
    if ~isfile(file_path)
        errordlg(['ファイルが見つかりません: ', file_path], 'ファイルエラー');
        return;
    end
    figure('WindowState', 'maximized');
    opts = detectImportOptions(file_path, 'NumHeaderLines', 2);
    T = readtable(file_path, opts);
    x = T{:,1};
    y = T{:,2};
    plot(x, y, 'LineWidth', 1.5);
    xlabel('Index');
    ylabel('Voltage [V]');
    title(sprintf('singleProbe/%s/%s', folder_date, file_name));
    grid on;
end

function answer = customDialog(default)
    d = dialog('Position', [500, 400, 350, 180], 'Name', 'singleProbe CSV選択');
    uicontrol('Parent', d, 'Style', 'text', 'Position', [20, 120, 100, 20], 'String', 'Date:');
    dateField = uicontrol('Parent', d, 'Style', 'edit', 'Position', [130, 120, 180, 25], 'String', default.dateValue);
    uicontrol('Parent', d, 'Style', 'text', 'Position', [20, 80, 100, 20], 'String', 'Shot:');
    shotField = uicontrol('Parent', d, 'Style', 'edit', 'Position', [130, 80, 180, 25], 'String', default.shotValue);
    btn = uicontrol('Parent', d, 'Position', [140, 25, 70, 30], 'String', 'OK', 'Callback', @btnCallback);
    answer = struct();
    function btnCallback(~, ~)
        answer.dateValue = get(dateField, 'String');
        answer.shotValue = get(shotField, 'String');
        delete(d);
    end
    uiwait(d);
    if isvalid(d)
        delete(d);
    end
end


if ~isdeployed && ~ismcc
    plot_singleProbe_main;
end
