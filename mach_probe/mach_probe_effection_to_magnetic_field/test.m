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
output_file = 'G:\\My Drive\\lab\\lab_data\\mach_probe\\mach_probe_effection_to_magnetic_field\\results.csv';

% Save the results to a CSV file
writetable(results_table, output_file);

% Display a message
fprintf('Results have been saved to %s\n', output_file);

%%
% CSVファイルのパスを指定
csv_file = 'G:\\My Drive\\lab\\lab_data\\mach_probe\\mach_probe_effection_to_magnetic_field\\results.csv';

% CSVファイルを読み込む
data_csv = readtable(csv_file);

% ファイル名から4桁の数値を抽出
extracted_numbers = cellfun(@(x) sscanf(x, 'G:\\My Drive\\lab\\lab_data\\pre_processed\\a039_%4d.mat'), data_csv.filename);

% Tからrlistとa039の数値を取得
date = 240610;
shotlist = [1:48];
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw'; % スプレッドシートのID
T = getTS6log(DOCID);
node = 'date';
T = searchlog(T, node, date);

% Extract valid a039 numbers and corresponding rlist values from T
a039_numbers = T.a039(shotlist);
valid_indices = ~isnan(a039_numbers);
a039_numbers_valid = a039_numbers(valid_indices);
rlist_valid = T.tripleProbeRPosition_cm_____________________0_5cm(shotlist(valid_indices)) * 0.01;

% Only consider data from data_csv that matches valid a039 numbers
[~, csv_indices] = ismember(a039_numbers_valid, extracted_numbers);

% Initialize figure
figure;

% データのプロット
plot(a039_numbers_valid, data_csv.psi_max(csv_indices), '-o', 'DisplayName', 'Psi Max');
hold on;
plot(a039_numbers_valid, data_csv.psi_min(csv_indices), '-o', 'DisplayName', 'Psi Min');
plot(a039_numbers_valid, data_csv.Bt_max(csv_indices), '-o', 'DisplayName', 'Bt Max');
plot(a039_numbers_valid, data_csv.Bt_min(csv_indices), '-o', 'DisplayName', 'Bt Min');

% Add shading based on rlist regions
colors = lines(length(rlist_valid)); % Use distinct colors for each region
for i = 1:length(rlist_valid)
    if i < length(rlist_valid)
        x_start = a039_numbers_valid(i);
        x_end = a039_numbers_valid(i+1);
        r_value = rlist_valid(i);

        % Apply shading for each region with different colors
        if r_value >= 0.05 && r_value < 0.1
            fill([x_start x_end x_end x_start], [ylim fliplr(ylim)], colors(i, :), ...
                'EdgeColor', colors(i, :), 'FaceAlpha', 0.3, 'DisplayName', sprintf('Region %.2f', r_value));
        end
    end
end

hold off;

% Set title and labels
title('Psi and Bt Max/Min Values with rlist Shading');
xlabel('File Number');
ylabel('Value');
legend('show');
grid on;
