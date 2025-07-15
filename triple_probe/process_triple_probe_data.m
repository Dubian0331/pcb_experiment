function [Te_matrix_eV, ne_matrix, R_values, time_values, matFileName] = process_triple_probe_data(date, filepath, DOCID, shotlist, params, case_name)
    % process_triple_probe_data: Te, neの計算と、ショット毎の生電流データをdouble行列で保存します。
    
    % ★★★ 追加1: 途中経過をプロットするための設定 ★★★
    enable_diagnostic_plots = false; % trueにするとプロットが有効になる
    plot_for_shot_list = [21, 33];    % ここにグラフを見たいショット番号を指定

    % --- パラメータ展開 ---
    V2 = params.V2; V3 = params.V3; A = params.A; S_probe = params.S_probe;
    time_us = params.time_us_unified;
    time_slice = params.time_slice_indices;
    lp_freq = params.lowpass_freq;
    smpl_freq = params.sampling_freq;
    sm_win = params.smoothing_window;
    final_sm_win = params.final_smoothing_window;

    % --- 物理定数 ---
    q = 1.60217663e-19; kb = 1.380649e-23; K2ev = 11604.5250061657; mi = A * 1.66054e-27;

    % --- スプレッドシートからR座標を取得 ---
    T = getTS6log(DOCID); T = searchlog(T, 'date', date);
    shot_indices = arrayfun(@(x) find(T.shot == x, 1), shotlist);
    r_list = T.tripleProbeRPosition_cm_(shot_indices) * 1e-2;

    % --- シンボリック計算 ---
    I1 = sym('I1'); I2 = sym('I2'); I3 = sym('I3'); Te = sym('Te'); ne_sym = sym('ne');
    eqn_Te = (I1+I2)/(I1+I3) == (1-exp(-q*V2/(kb*Te)))/(1-exp(-q*V3/(kb*Te)));
    Te_symbolic = solve(eqn_Te, Te);
    eqn_ne = exp(-0.5)*S_probe*ne_sym*q*sqrt(kb*Te/mi) == (I3-I2*exp(-q*(V3-V2)/(kb*Te)))/(1-exp(-q*(V3-V2)/(kb*Te)));
    ne_symbolic = solve(eqn_ne, ne_sym);

    % --- ショットごとの処理 ---
    num_shots = length(shotlist);
    Te_matrix_shots = zeros(num_shots, length(time_us));
    ne_matrix_shots = zeros(num_shots, length(time_us));
    valid_shots_mask = false(1, num_shots);
    
    % ★★★ 変更点1: double型の行列をゼロで初期化 ★★★
    num_time_points = length(time_slice);
    time_raw_all_shots = zeros(num_time_points, num_shots);
    I2_raw_all_shots = zeros(num_time_points, num_shots);
    I3_raw_all_shots = zeros(num_time_points, num_shots);

    for idx = 1:num_shots
        shot = shotlist(idx);
        filename = fullfile(filepath, num2str(date), ['ES_', num2str(date), sprintf('%03d', shot), '.csv']);
        if ~exist(filename, 'file'), continue; end
        
        valid_shots_mask(idx) = true;
        data = readmatrix(filename);
        
        time_raw = data(time_slice, 1);
        I2_raw = data(time_slice, 36);
        I3_raw = data(time_slice, 37);
        
        % ★★★ 変更点2: 行列の対応する列にデータを格納 ★★★
        time_raw_all_shots(:, idx) = time_raw;
        I2_raw_all_shots(:, idx) = I2_raw;
        I3_raw_all_shots(:, idx) = I3_raw;

        % --- Te/ne計算のためのデータ前処理 ---
        % 負の値をNaNに置換
        I2_raw(I2_raw < 0) = NaN;
        I3_raw(I3_raw < 0) = NaN;
        
        % NaNを線形補間し、lowpassフィルタを適用
        I2_processed = lowpass(fillmissing(I2_raw, 'linear'), lp_freq, smpl_freq);
        I3_processed = lowpass(fillmissing(I3_raw, 'linear'),lp_freq, smpl_freq);
        
        % スムージング
        I2_processed = smoothdata(I2_processed, 'movmedian', sm_win);
        I3_processed = smoothdata(I3_processed, 'movmedian', sm_win);

        % フィルタ後のデータに、元のNaNの位置を復元
        I2_processed(isnan(I2_raw)) = NaN;
        I3_processed(isnan(I3_raw)) = NaN;
    
        % 信号が小さすぎる区間を無効化するための閾値処理
        peak_current = max(max(I2_processed, [], 'omitnan'), max(I3_processed, [], 'omitnan'));
        if ~isempty(peak_current) && peak_current > 0
            current_threshold = peak_current * 0.05; % ピーク値の5%を閾値とする
            I2_processed(I2_processed < current_threshold) = NaN;
            I3_processed(I3_processed < current_threshold) = NaN;
        end

        I1_processed = I2_processed + I3_processed;

        try
            Te_values_K = calculate_Te(Te_symbolic, I1_processed, I2_processed, I3_processed, K2ev);
        catch
            disp(['An error occurred during Te calculation for shot ', num2str(shot), ':']);
            disp(ME.message);
            Te_values_K = NaN(size(I1_processed));
        end
        try
            ne_values = real(double(subs(ne_symbolic, {Te, I2, I3}, {Te_values_K, I2_processed, I3_processed})));
            ne_values(ne_values < 0) = NaN;
        catch
            ne_values = NaN(size(I2_processed));
        end
        
        % --- ★★★ ここからが修正・追加箇所 ★★★ ---
        if enable_diagnostic_plots && ismember(shot, plot_for_shot_list)
            % Figure 1: 電流値のプロット
            figure('Name', ['Diagnostic Plot for Shot ', num2str(shot)]);
            hold on;
            grid on;
            plot(time_raw, I2_raw, ':', 'Color', [1 0.6 0.6], 'DisplayName', 'I2 (Raw)');
            plot(time_raw, I3_raw, ':', 'Color', [0.6 0.6 1], 'DisplayName', 'I3 (Raw)');
            plot(time_raw, I2_processed, '-', 'Color', 'r', 'LineWidth', 1.5, 'DisplayName', 'I2 (Processed)');
            plot(time_raw, I3_processed, '-', 'Color', 'b', 'LineWidth', 1.5, 'DisplayName', 'I3 (Processed)');
            hold off;
            title(['Current Processing for Shot ', num2str(shot)]);
            xlabel('Time (us)');
            ylabel('Current [A]');
            xlim([min(time_raw), max(time_raw)]);
            legend;

            % Figure 2: 電流比とTeのプロット
            figure('Name', ['Te & Ratio for Shot ', num2str(shot)]);
            
            % 左側のY軸に電流比をプロット
            yyaxis left
            Te_calc_ratio = (I1_processed + I2_processed) ./ (I1_processed + I3_processed);
            plot(time_raw, Te_calc_ratio, 'k-');
            ylabel('Current Ratio (I1+I2)/(I1+I3)');
            ylim([0, 2]); % 電流比の表示範囲
            
            hold on;
            grid on;
            
            % 右側のY軸にTe [eV]をプロット
            yyaxis right
            plot(time_raw, Te_values_K / K2ev, 'm-', 'LineWidth', 1.5);
            ylabel('Electron Temperature [eV]');
            ylim([0, 25]); % Teの表示範囲
            
            hold off;
            
            title(['Calculated Te and Current Ratio for Shot ', num2str(shot)]);
            xlabel('Time (us)');
            xlim([min(time_raw), max(time_raw)]);
            legend({'Current Ratio', 'Te [eV]'}, 'Location', 'northeast');

            drawnow;
        end

        Te_interp = interp1(time_raw, Te_values_K, time_us, 'linear', NaN);
        ne_interp = interp1(time_raw, ne_values, time_us, 'linear', NaN);
        Te_matrix_shots(idx, :) = smoothdata(Te_interp, 'movmedian', final_sm_win);
        ne_matrix_shots(idx, :) = smoothdata(ne_interp, 'movmedian', final_sm_win);
    end

    % --- 集計処理 ---
    % 有効なショットのデータのみを抽出
    Te_matrix = Te_matrix_shots(valid_shots_mask, :);
    ne_matrix = ne_matrix_shots(valid_shots_mask, :);
    r_list_valid = r_list(valid_shots_mask);
    
    % ★★★ 変更点3: 有効なショットに対応する列のみを抽出 ★★★
    time_raw_per_shot = time_raw_all_shots(:, valid_shots_mask);
    I2_raw_per_shot = I2_raw_all_shots(:, valid_shots_mask);
    I3_raw_per_shot = I3_raw_all_shots(:, valid_shots_mask);

    for i = 1:size(Te_matrix, 1)
        Te_matrix(i, :) = fillmissing(Te_matrix(i, :), 'linear', 'EndValues', 'none');
        ne_matrix(i, :) = fillmissing(ne_matrix(i, :), 'linear', 'EndValues', 'none');
    end

    % 同じr_listの値でTeとneの平均値を計算
    R_values= unique(r_list_valid); % r_list_validを使用
    Te_avg_matrix = zeros(length(R_values), length(time_us));
    ne_avg_matrix = zeros(length(R_values), length(time_us));
    for i = 1:length(R_values)
        idx_r = (r_list_valid == R_values(i)); % r_list_validを使用
        Te_avg_matrix(i, :) = mean(Te_matrix(idx_r, :), 1, 'omitnan');
        ne_avg_matrix(i, :) = mean(ne_matrix(idx_r, :), 1, 'omitnan');
    end
    
    % 平均化後の時間方向の最終スムージング
    for i = 1:size(Te_avg_matrix, 1)
        Te_avg_matrix(i, :) = smoothdata(Te_avg_matrix(i, :), 'movmedian', 7, 'omitnan'); 
        ne_avg_matrix(i, :) = smoothdata(ne_avg_matrix(i, :), 'movmedian', 7, 'omitnan');
    end
    
    % 各時刻のデータ点について、R方向（列方向）のNaNを線形補間します。
    Te_avg_matrix = fillmissing(Te_avg_matrix, 'linear', 1, 'EndValues', 'none');
    ne_avg_matrix = fillmissing(ne_avg_matrix, 'linear', 1, 'EndValues', 'none');
    
    % さらに、R方向に平滑化をかけ、より滑らかな分布にします。
    Te_avg_matrix = smoothdata(Te_avg_matrix, 1, 'movmedian', 3, 'omitnan');
    ne_avg_matrix = smoothdata(ne_avg_matrix, 1, 'movmedian', 3, 'omitnan');


    % for i = 1:size(Te_matrix, 1)
    %     Te_matrix(i, :) = fillmissing(Te_matrix(i, :), 'linear', 'EndValues', 'nearest'); % 端のNaNも埋める
    %     ne_matrix(i, :) = fillmissing(ne_matrix(i, :), 'linear', 'EndValues', 'nearest');
    % end
    % 
    % R_values = unique(r_list_valid);
    % Te_avg_matrix = zeros(length(R_values), length(time_us));
    % ne_avg_matrix = zeros(length(R_values), length(time_us));
    % for i = 1:length(R_values)
    %     idx_r = (r_list_valid == R_values(i));
    %     Te_avg_matrix(i, :) = mean(Te_matrix(idx_r, :), 1, 'omitnan');
    %     ne_avg_matrix(i, :) = mean(ne_matrix(idx_r, :), 1, 'omitnan');
    % end
    % 
    % for i = 1:size(Te_avg_matrix, 1) % 時間方向
    %     Te_avg_matrix(i, :) = smoothdata(Te_avg_matrix(i, :), 'movmedian', 7, 'omitnan'); 
    %     ne_avg_matrix(i, :) = smoothdata(ne_avg_matrix(i, :), 'movmedian', 7, 'omitnan');
    % end
    
    Te_matrix_eV = Te_avg_matrix / K2ev;
    ne_matrix = ne_avg_matrix;
    time_values = time_us;
    valid_shotlist = shotlist(valid_shots_mask);

    % --- MATファイル保存 ---
    targetFolderPath = fullfile(filepath, '..', 'triple_probe', 'mat', num2str(date));
    if ~exist(targetFolderPath, 'dir'), mkdir(targetFolderPath); end
    matFileName = fullfile(targetFolderPath, ['triple_data2D_', case_name, '.mat']);
    
    triple_data2D = struct();
    triple_data2D.Te_avg = reshape(Te_matrix_eV, [length(R_values), 1, length(time_values)]);
    triple_data2D.ne_avg = reshape(ne_matrix, [length(R_values), 1, length(time_values)]);
    
    % ★★★ 変更点4: 保存する変数リストを更新 ★★★
    save(matFileName, ...
        'triple_data2D', ...
        'R_values', ...
        'time_values', ...
        'shotlist', ...
        'valid_shotlist', ...
        'time_raw_per_shot', ...  % ショット毎の生の時間データ (double行列)
        'I2_raw_per_shot', ...    % ショット毎の生のI2データ (double行列)
        'I3_raw_per_shot', ...    % ショット毎の生のI3データ (double行列)
        'time_slice');

    disp(['MATファイル ', matFileName, ' を保存しました。']);
end

function Te_values = calculate_Te(Te_symbolic, I1, I2, I3, K2ev)
    I1_c = I1; I2_c = I2; I3_c = I3;
    ratio_arg = (I1_c + I3_c) ./ (I2_c - I3_c);
    problem_idx = abs(I2_c - I3_c) < 1e-12 | ratio_arg <= 0 | abs(ratio_arg - 1) < 1e-12;
    I1_c(problem_idx) = NaN;
    I2_c(problem_idx) = NaN;
    I3_c(problem_idx) = NaN;
    Te_values = real(double(subs(Te_symbolic, {'I1', 'I2', 'I3'}, {I1_c, I2_c, I3_c})));
    Te_values(Te_values < 0 | Te_values > 100 * K2ev) = NaN;
end