% Define the file name pattern
file_pattern = 'G:\\My Drive\\lab\\lab_data\\pre_processed\\a039_4%d.mat';

% Define the range of file numbers
file_numbers = 289:319;

% Initialize an empty cell array to store the results
results = {'filename', 'psi_max', 'psi_min', 'Bt_max', 'Bt_min'};

% Loop through each file number
for i = file_numbers
    % Create the full file path
    file_name = sprintf(file_pattern, i);
    
    % Check if the file exists
    if exist(file_name, 'file')
        % Load the .mat file
        data = load(file_name);
        
        % Extract 'psi' and 'Bt' from the loaded data structure
        % Replace 'data2D' with the actual variable name in the .mat file
        if isfield(data, 'data2D')
            psi_data = data.data2D.psi; % Adjust if needed
            Bt_data = data.data2D.Bt;   % Adjust if needed

            % Calculate the max and min for psi and Bt
            psi_max = max(psi_data(:));
            psi_min = min(psi_data(:));
            Bt_max = max(Bt_data(:));
            Bt_min = min(Bt_data(:));

            % Append the results to the cell array
            results(end+1, :) = {file_name, psi_max, psi_min, Bt_max, Bt_min};
        else
            fprintf('Expected field "data2D" not found in %s. Skipping...\n', file_name);
            continue;
        end
    else
        % Display a warning and continue if the file does not exist
        fprintf('File %s does not exist. Skipping...\n', file_name);
        continue;
    end
end

% Convert the cell array to a table
results_table = cell2table(results(2:end, :), 'VariableNames', results(1, :));

% Define the output file path
output_file = 'G:\My Drive\lab\lab_data\mach_probe\mach_probe_effection_to_magnetic_field\results.csv';

% Save the results to a CSV file
writetable(results_table, output_file);

% Display a message
fprintf('Results have been saved to %s\n', output_file);


%%
% CSVファイルのパスを指定
csv_file = 'G:\\My Drive\\lab\\lab_data\\mach_probe\\mach_probe_effection_to_magnetic_field\\results.csv';

% CSVファイルを読み込む
data = readtable(csv_file);

% ファイル名から4桁の数値を抽出
extracted_numbers = cellfun(@(x) sscanf(x, 'G:\\My Drive\\lab\\lab_data\\pre_processed\\a039_%4d.mat'), data.filename);

% フィギュアを作成
figure;

% 1番目のサブプロット: psiの最大値
subplot(2,2,1);
plot(extracted_numbers, data.psi_max, '-o');
title('Psi Max');
xlabel('File Number');
ylabel('Psi Max');
grid on;

% 2番目のサブプロット: psiの最小値
subplot(2,2,2);
plot(extracted_numbers, data.psi_min, '-o');
title('Psi Min');
xlabel('File Number');
ylabel('Psi Min');
grid on;

% 3番目のサブプロット: Btの最大値
subplot(2,2,3);
plot(extracted_numbers, data.Bt_max, '-o');
title('Bt Max');
xlabel('File Number');
ylabel('Bt Max');
grid on;

% 4番目のサブプロット: Btの最小値
subplot(2,2,4);
plot(extracted_numbers, data.Bt_min, '-o');
title('Bt Min');
xlabel('File Number');
ylabel('Bt Min');
grid on;

% フィギュア全体のタイトル
sgtitle('Psi and Bt Max/Min Values by File Number');

%%
% CSVファイルのパスを指定
csv_file = 'G:\\My Drive\\lab\\lab_data\\mach_probe\\mach_probe_effection_to_magnetic_field\\results.csv';

% CSVファイルを読み込む
data = readtable(csv_file);

% ファイル名から4桁の数値を抽出
extracted_numbers = cellfun(@(x) sscanf(x, 'G:\\My Drive\\lab\\lab_data\\pre_processed\\a039_%4d.mat'), data.filename);

% フィギュアを作成
figure;

% Define shading regions (adjust the range to fit your data)
shading_regions = [
    struct('x_start', 4290, 'x_end', 4300, 'color', [0.9 0.9 0.9]);
    struct('x_start', 4305, 'x_end', 4315, 'color', [0.8 0.8 0.8])
];

% 1番目のサブプロット: psiの最大値
subplot(2,2,1);
plot(extracted_numbers, data.psi_max, '-o');
title('Psi Max');
xlabel('File Number');
ylabel('Psi Max');
grid on;
% Get the current y-limits
y_limits = ylim;
% Add common shading
hold on;
for i = 1:length(shading_regions)
    fill([shading_regions(i).x_start shading_regions(i).x_end shading_regions(i).x_end shading_regions(i).x_start], ...
         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], shading_regions(i).color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
hold off;

% 2番目のサブプロット: psiの最小値
subplot(2,2,2);
plot(extracted_numbers, data.psi_min, '-o');
title('Psi Min');
xlabel('File Number');
ylabel('Psi Min');
grid on;
% Get the current y-limits
y_limits = ylim;
% Add common shading
hold on;
for i = 1:length(shading_regions)
    fill([shading_regions(i).x_start shading_regions(i).x_end shading_regions(i).x_end shading_regions(i).x_start], ...
         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], shading_regions(i).color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
hold off;

% 3番目のサブプロット: Btの最大値
subplot(2,2,3);
plot(extracted_numbers, data.Bt_max, '-o');
title('Bt Max');
xlabel('File Number');
ylabel('Bt Max');
grid on;
% Get the current y-limits
y_limits = ylim;
% Add common shading
hold on;
for i = 1:length(shading_regions)
    fill([shading_regions(i).x_start shading_regions(i).x_end shading_regions(i).x_end shading_regions(i).x_start], ...
         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], shading_regions(i).color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
hold off;

% 4番目のサブプロット: Btの最小値
subplot(2,2,4);
plot(extracted_numbers, data.Bt_min, '-o');
title('Bt Min');
xlabel('File Number');
ylabel('Bt Min');
grid on;
% Get the current y-limits
y_limits = ylim;
% Add common shading
hold on;
for i = 1:length(shading_regions)
    fill([shading_regions(i).x_start shading_regions(i).x_end shading_regions(i).x_end shading_regions(i).x_start], ...
         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], shading_regions(i).color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
hold off;

% フィギュア全体のタイトル
sgtitle('Psi and Bt Max/Min Values by File Number with Shading');


%%
% 提供されたデータ（ここでは省略）
data = simple_list;

% シェーディング用の構造体を作成
shading_regions = [];
i = find(data(:, 1) == 4280);  % 4289に対応するインデックスを取得
while i < size(data, 1)
    x_start = data(i, 1);
    color_value = data(i, 2);
    
    % 同じ値が続く区間を探す
    j = i;
    while j < size(data, 1) && data(j + 1, 2) == color_value
        j = j + 1;
    end
    
    x_end = data(j, 1);
    
    % カラーを設定（グレーと白）
    if mod(i, 2) == 1
        color = [0.5, 0.5, 0.5]; % グレー
    else
        color = [0, 0, 0]; % 白
    end
    
    shading_regions = [shading_regions; struct('x_start', x_start, 'x_end', x_end, 'color', color)];
    
    % インデックスを次の区間のスタートに移動
    i = j + 1;
end

% CSVファイルのパスを指定
csv_file = 'G:\\My Drive\\lab\\lab_data\\mach_probe\\mach_probe_effection_to_magnetic_field\\results.csv';

% CSVファイルを読み込む
data_csv = readtable(csv_file);

% ファイル名から4桁の数値を抽出
extracted_numbers = cellfun(@(x) sscanf(x, 'G:\\My Drive\\lab\\lab_data\\pre_processed\\a039_%4d.mat'), data_csv.filename);

% フィギュアを作成
figure;

% 1番目のサブプロット: psiの最大値
subplot(2,2,1);
plot(extracted_numbers, data_csv.psi_max, '-o');
title('Psi Max');
xlabel('File Number');
ylabel('Psi Max');
grid on;
y_limits = ylim;
hold on;
for k = 1:length(shading_regions)
    fill([shading_regions(k).x_start shading_regions(k).x_end shading_regions(k).x_end shading_regions(k).x_start], ...
         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], shading_regions(k).color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
hold off;

% 2番目のサブプロット: psiの最小値
subplot(2,2,2);
plot(extracted_numbers, data_csv.psi_min, '-o');
title('Psi Min');
xlabel('File Number');
ylabel('Psi Min');
grid on;
y_limits = ylim;
hold on;
for k = 1:length(shading_regions)
    fill([shading_regions(k).x_start shading_regions(k).x_end shading_regions(k).x_end shading_regions(k).x_start], ...
         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], shading_regions(k).color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
hold off;

% 3番目のサブプロット: Btの最大値
subplot(2,2,3);
plot(extracted_numbers, data_csv.Bt_max, '-o');
title('Bt Max');
xlabel('File Number');
ylabel('Bt Max');
grid on;
y_limits = ylim;
hold on;
for k = 1:length(shading_regions)
    fill([shading_regions(k).x_start shading_regions(k).x_end shading_regions(k).x_end shading_regions(k).x_start], ...
         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], shading_regions(k).color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
hold off;

% 4番目のサブプロット: Btの最小値
subplot(2,2,4);
plot(extracted_numbers, data_csv.Bt_min, '-o');
title('Bt Min');
xlabel('File Number');
ylabel('Bt Min');
grid on;
y_limits = ylim;
hold on;
for k = 1:length(shading_regions)
    fill([shading_regions(k).x_start shading_regions(k).x_end shading_regions(k).x_end shading_regions(k).x_start], ...
         [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], shading_regions(k).color, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
end
hold off;

% フィギュア全体のタイトル
sgtitle('Psi and Bt Max/Min Values by File Number with Shading');
saveas(gcf, '')
