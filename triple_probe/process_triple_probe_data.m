function [Te_matrix_eV, ne_matrix, R_values, time_values, matFileName] = process_triple_probe_data(date, filepath, DOCID, shotlist, params, case_name)
    % process_triple_probe_data: 元のスクリプトのロジックを忠実に再現し、Teとneを計算します。
    
    % --- パラメータ展開 ---
    V2 = params.V2; V3 = params.V3; A = params.A; S_probe = params.S_probe;
    time_us = params.time_us_unified;
    time_slice = params.time_slice_indices;
    noise_thr = params.current_noise_threshold;
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
    I1 = sym('I1');
    I2 = sym('I2');
    I3 = sym('I3');
    Te = sym('Te');
    ne_sym = sym('ne'); % 'ne'から'ne_sym'へ変更

    eqn_Te = (I1+I2)/(I1+I3) == (1-exp(-q*V2/(kb*Te)))/(1-exp(-q*V3/(kb*Te)));
    Te_symbolic = solve(eqn_Te, Te);

    eqn_ne = exp(-0.5)*S_probe*ne_sym*q*sqrt(kb*Te/mi) == (I3-I2*exp(-q*(V3-V2)/(kb*Te)))/(1-exp(-q*(V3-V2)/(kb*Te)));
    ne_symbolic = solve(eqn_ne, ne_sym);

    % --- ショットごとの処理 ---
    num_shots = length(shotlist);
    Te_matrix_shots = zeros(num_shots, length(time_us));
    ne_matrix_shots = zeros(num_shots, length(time_us));
    valid_shots_mask = false(1, num_shots);

    for idx = 1:num_shots
        shot = shotlist(idx);
        filename = fullfile(filepath, num2str(date), ['ES_', num2str(date), sprintf('%03d', shot), '.csv']);
        if ~exist(filename, 'file'), continue; end
        
        valid_shots_mask(idx) = true;
        data = readmatrix(filename);
        time = data(time_slice, 1);
        I2_values = data(time_slice, 36);
        I3_values = data(time_slice, 37);

        I2_values(I2_values < 0 | abs(I2_values) < noise_thr) = NaN;
        I3_values(I3_values < 0 | abs(I3_values) < noise_thr) = NaN;

        I2_values = fillmissing(I2_values, 'linear');
        I3_values = fillmissing(I3_values, 'linear');
        
        I2_values = lowpass(I2_values, lp_freq, smpl_freq);
        I3_values = lowpass(I3_values, lp_freq, smpl_freq);
        
        I2_values = smoothdata(I2_values, 'movmedian', sm_win);
        I3_values = smoothdata(I3_values, 'movmedian', sm_win);
        I1_values = I2_values + I3_values;

        % Te 計算 (元のロジックを流用)
        try
            Te_values_K = calculate_Te(Te_symbolic, I1_values, I2_values, I3_values, K2ev);
        catch
            Te_values_K = NaN(size(I1_values));
        end

        % ne 計算
        try
            ne_values = real(double(subs(ne_symbolic, {Te, I2, I3}, {Te_values_K, I2_values, I3_values})));
            ne_values(ne_values < 0) = NaN;
        catch
            ne_values = NaN(size(I2_values));
        end
        
        % 時間軸統一と最終平滑化
        Te_matrix_shots(idx, :) = interp1(time, Te_values_K, time_us, 'linear', 0);
        ne_matrix_shots(idx, :) = interp1(time, ne_values, time_us, 'linear', 0);
        
        Te_matrix_shots(idx, :) = smoothdata(Te_matrix_shots(idx, :), 'movmedian', final_sm_win);
        ne_matrix_shots(idx, :) = smoothdata(ne_matrix_shots(idx, :), 'movmedian', final_sm_win);
    end

    % --- 集計処理 ---
    Te_matrix = Te_matrix_shots(valid_shots_mask, :);
    ne_matrix = ne_matrix_shots(valid_shots_mask, :);
    r_list_valid = r_list(valid_shots_mask);

    for i = 1:size(Te_matrix, 1) % 元のスクリプトのNaN埋め処理
        Te_matrix(i, :) = fillmissing(Te_matrix(i, :), 'linear');
        ne_matrix(i, :) = fillmissing(ne_matrix(i, :), 'linear');
    end

    R_values = unique(r_list_valid);
    Te_avg_matrix = zeros(length(R_values), length(time_us));
    ne_avg_matrix = zeros(length(R_values), length(time_us));

    for i = 1:length(R_values)
        idx_r = (r_list_valid == R_values(i));
        Te_avg_matrix(i, :) = mean(Te_matrix(idx_r, :), 1, 'omitnan');
        ne_avg_matrix(i, :) = mean(ne_matrix(idx_r, :), 1, 'omitnan');
    end
    
    Te_matrix_eV = Te_avg_matrix / K2ev;
    ne_matrix = ne_avg_matrix;
    time_values = time_us;

    % --- MATファイル保存 ---
    targetFolderPath = fullfile(filepath, '..', 'triple_probe', 'mat', num2str(date));
    if ~exist(targetFolderPath, 'dir'), mkdir(targetFolderPath); end
    matFileName = fullfile(targetFolderPath, ['triple_data2D_', case_name, '.mat']);
    triple_data2D.Te = reshape(Te_matrix_eV, [length(R_values), 1, length(time_values)]);
    triple_data2D.ne = reshape(ne_matrix, [length(R_values), 1, length(time_values)]);
    save(matFileName, 'triple_data2D', 'R_values', 'time_values');
end

function Te_values = calculate_Te(Te_symbolic, I1, I2, I3, K2ev)
    % Te計算のロバストなサブルーチン
    I1_c = I1; I2_c = I2; I3_c = I3;
    ratio_arg = (I1_c + I3_c) ./ (I2_c - I3_c);
    problem_idx = abs(I2_c - I3_c) < 1e-12 | ratio_arg <= 0 | abs(ratio_arg - 1) < 1e-12;
    I1_c(problem_idx) = NaN;
    I2_c(problem_idx) = NaN;
    I3_c(problem_idx) = NaN;
    Te_values = real(double(subs(Te_symbolic, {'I1', 'I2', 'I3'}, {I1_c, I2_c, I3_c})));
    Te_values(Te_values < 0 | Te_values > 100 * K2ev) = NaN;
end