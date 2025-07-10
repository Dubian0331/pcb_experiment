%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Et profile with error bars (Time Evolution)
% Multiple Files with Different Sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% ショットリスト
cases = {struct('name', 'CaseI', 'shots', 4289:4319, 'color', 'r', 'col', 22), ...
         struct('name', 'CaseO', 'shots', 4363:4389, 'color', 'b', 'col', 20)};

% 時間の範囲を設定 (例: 470 µs から 600 µs)
timeRangeStart = 470; % 開始時間
timeRangeEnd = 600;   % 終了時間

% グラフ用の初期化
figure;
hold on;
grid on;
xlabel('Time [\mus]');
ylabel('E_t [V/m]');
% xlim([460 510])
yline(0, 'k--', 'HandleVisibility', 'off');
title('Comparison of Case-I and Case-O');

% 各ケースのデータを計算
for i = 1:length(cases)
    caseData = cases{i};
    list_shots = caseData.shots;
    Et_timeSeries_all = [];
    time_all = [];

    for shot = list_shots
        % データをロード
        fileName = sprintf('G:\\My Drive\\lab\\lab_data\\pre_processed\\a039_%d.mat', shot);
        try
            load(fileName);

            % データ構造の確認
            if ~isfield(data2D, 'psi') || ~isfield(data2D, 'trange')
                fprintf('Skipping shot %d: data2D structure is invalid.\n', shot);
                continue;
            end

            % 配列サイズを確認
            [rows, cols, depth] = size(data2D.psi);
            if ~(rows == 20 && cols == 40 && depth == 201)
                fprintf('Skipping shot %d: data2D.psi size is not 20x40x201.\n', shot);
                continue;
            end

            % 時間インデックスの抽出
            timeIndices = find(data2D.trange >= timeRangeStart & data2D.trange <= timeRangeEnd);
            time = data2D.trange(timeIndices);

            % 初回のみ時間を保存
            if isempty(time_all)
                time_all = time(1:end-1); % `diff` で長さが減るため調整
            end

            % 指定された列データを取得
            col = caseData.col; % CaseI: 22列目, CaseO: 20列目
            Et_timeSeries = zeros(1, length(timeIndices) - 1);
            for idx = 1:length(timeIndices) - 1
                t1 = timeIndices(idx);
                t2 = timeIndices(idx + 1);
                
                % psi の差分から Et を計算
                psi_diff = data2D.psi(:, col, t2) - data2D.psi(:, col, t1);
                dt = data2D.trange(t2) - data2D.trange(t1);
                Et_slice = -psi_diff / (2 * pi * grid2D.rq(col)); % Et 計算

                % 平均値を保存
                Et_timeSeries(idx) = mean(Et_slice); % 平均値を計算
            end

            % 計算結果を保存
            Et_timeSeries_all = [Et_timeSeries_all; Et_timeSeries];
        catch ME
            fprintf('Error processing shot %d: %s\n', shot, ME.message);
            continue;
        end
    end

    % 平均と標準偏差の計算
    Et_mean = mean(Et_timeSeries_all, 1);
    Et_std = std(Et_timeSeries_all, 0, 1);

    % エラーバー付きプロット
    errorbar(time_all, Et_mean, Et_std, caseData.color, 'LineWidth', 1.5, 'DisplayName', caseData.name);
end

legend;
hold off;

% 結果の保存
saveDir = 'G:\\My Drive\\lab\\lab_data\\mach_probe\\Et_profile\\';
fileName = 'Comparison_Case-I_Case-O_Et_time_evolution';
format = 'png';
fig = gcf;
saveas(fig, fullfile(saveDir, [fileName '.' format]));
