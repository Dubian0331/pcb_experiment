%%%%%%%%%%%%%%%%%%%%%%%%%%%
% memo
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Case-I:6/10, Xpoint:-0.01m -> val(;, 20, target_time) & val(:, 7, target_time), Target time: "400 + target_time". 499us
% startFile = 4289; 
% endFile = 4319;   
% X_point_is_here = 0.24;
% target_time = 99;

% Case-O:6/11, Xpoint:0.03m -> val(;, 22, target_time) & val(:, 9, target_time), file_start:4363, file_end:4389
% startFile = 4363; % 開始ファイル番号
% endFile = 4389;   % 終了ファイル番号
% X_point_is_here = 0.20;
% target_time = 90;
% target_time = 112;

%%
% --- プロットするケースを選択 ---
% Case-I
% date = 240610;
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240610018.mat');
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240610\grad_P_data_240610.mat');
% z_index = 20; % プロットするzの位置 (z=0に対応)

% % case-O
date = 240611;
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240611055.mat');
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240611\grad_P_data_240611.mat');
z_index = 22; % プロットするzの位置 (z=0に対応)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time_us = 480; % プロットを開始する時刻 [us]
time_step_us = 1;    % プロットの時間間隔 [us]
num_plots_y = 5;     % 縦に並べるプロット数
num_plots_x = 6;     % 横に並べるプロット数
dispGradP = false; % gradPを表示するか
dispSum = false; % J×BとgradPの合計値を表示するか
isSave = false; % 画像を保存するかどうか
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- X点の位置の取得 ---
% get_axis_x_multi.m がパス上にある必要があります
fprintf('Calculating X-point locations for all time steps...\n');
[~, xPointList] = get_axis_x_multi(grid2D, data2D);
fprintf('X-point calculation complete.\n');

%% --- プロット処理 ---
figure('WindowState', 'maximized');
r_axis_jxb = grid2D.rq(:,1);
jxb_time_axis = data2D.trange;
total_plots = num_plots_y * num_plots_x;

for plot_idx = 1:total_plots
    time_in_us = start_time_us + (plot_idx - 1) * time_step_us;
    
    subplot(num_plots_y, num_plots_x, plot_idx);
    hold on;
    grid on;
    
    % --- JxBの計算 ---
    [~, target_time] = min(abs(jxb_time_axis - time_in_us));
    Jt = data2D.Jt(:, z_index, target_time);
    Bz = data2D.Bz(:, z_index, target_time);
    jxb_to_plot = (Jt .* Bz) / 1000;
    
    plot(r_axis_jxb, jxb_to_plot, 'b-', 'LineWidth', 1, 'DisplayName', 'J_tB_z');
    
    % --- grad Pの計算とプロット ---
    if dispGradP || dispSum
        [~, grad_P_time_idx] = min(abs(time_axis - time_in_us));
        grad_P_to_plot = grad_P_2D(grad_P_time_idx, :);
        grad_P_aligned = interp1(R_axis, -grad_P_to_plot, r_axis_jxb, 'pchip', 'extrap') / 1000;
        
        if dispGradP
            plot(r_axis_jxb, grad_P_aligned, 'r', 'LineWidth', 1, 'DisplayName', '-grad P');
        end
        if dispSum
            force_sum = jxb_to_plot + grad_P_aligned;
            plot(r_axis_jxb, force_sum, 'k--', 'LineWidth', 1, 'DisplayName', 'Sum');
        end
    end

    % --- X点について ---
    x_point_r = xPointList.r(target_time);
    
    % X点が見つかっている場合（NaNでない場合）のみプロット
    if ~isnan(x_point_r)
        % 凡例が重複しないように工夫
        if plot_idx == 1
            xline(x_point_r, 'k:', 'LineWidth', 1.5, 'DisplayName', 'X-point');
        else
            xline(x_point_r, 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
    end
    
    % --- グラフ装飾 ---
    yline(0, 'k-', 'HandleVisibility', 'off');
    title(['t = ', num2str(time_in_us), ' us']);
    if plot_idx == 1, legend('FontSize', 7); end
    if plot_idx > total_plots - num_plots_x, xlabel('r [m]'); end
    if mod(plot_idx - 1, num_plots_x) == 0
        ylabel('Force [kN/m^3]');
    end
end
sgtitle(['Force Balance Analysis for Date: ', num2str(date)], 'FontSize', 16, 'FontWeight', 'bold');

% --- 画像保存処理 ---
if isSave
    % 保存先ディレクトリ
    base_dir = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\j×B_calc\\figure';
    date_dir = fullfile(base_dir, num2str(date));
    if ~exist(date_dir, 'dir')
        mkdir(date_dir);
    end

    % ファイル名作成
    if date == 240610
        case_name = 'CaseI';
    elseif date == 240611
        case_name = 'CaseO';
    else
        error('Unsupported date. Please set date to 240610 or 240611.');
    end

    start_str = num2str(start_time_us);
    end_str = num2str(start_time_us + (total_plots-1)*time_step_us);
    save_name = sprintf('ForceBalance_%s_%s_%sus-%sus.png', case_name, num2str(date), start_str, end_str);
    save_path = fullfile(date_dir, save_name);

    saveas(gcf, save_path);
    disp(['Figure saved to: ', save_path]);
end