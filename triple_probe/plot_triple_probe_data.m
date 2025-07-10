function plot_triple_probe_data(Te_data_2D, ne_data_2D, R_values, time_values, options)
% plot_triple_probe_data: Teとneの2次元コンタープロットを作成します。
    
    [T_grid, R_grid] = meshgrid(time_values, R_values);
    
    % Teプロット
    fig_Te = figure('Name', [options.case_name, ' Te Plot']);
    contourf(T_grid, R_grid, Te_data_2D, 100, 'LineColor', 'none');
    cb = colorbar;
    cb.Title.String = 'T_{e} [eV]';
    cb.Title.FontSize = 10; 
    cb.FontSize = 10;  
    clim(options.Te_clim); xlim(options.time_lim); ylim(options.R_lim);
    colormap('jet');

    title([options.case_name, ' plot of T_{e} [eV]']); xlabel('Time [\mus]'); ylabel('R [m]');
    set(gca, 'FontSize', 12);

    % neプロット
    fig_ne = figure('Name', [options.case_name, ' ne Plot']);
    contourf(T_grid, R_grid, ne_data_2D, 100, 'LineColor', 'none');
    cb = colorbar;
    cb.Title.String = 'n_{e} [m^{-3}]';
    cb.Title.FontSize = 10; 
    cb.FontSize = 10;  
    clim(options.ne_clim); xlim(options.time_lim); ylim(options.R_lim);
    
    colormap('jet');
    title([options.case_name, ' plot of n_{e} [m^{-3}]']); xlabel('Time [\mus]'); ylabel('R [m]');
    set(gca, 'FontSize', 12);

    % --- グラフの保存 ---
    if options.save_figures
        saveDir = options.save_path;
        if ~exist(saveDir, 'dir'), mkdir(saveDir); end
        
        % ファイル名に日付とケース名を含める  例: "240611_case-o"
        base_filename = [num2str(options.date), '_', lower(strrep(options.case_name, ' ', '_'))];

        % Teグラフのフルパスを生成して保存
        filename_Te = fullfile(saveDir, [base_filename, '_te.png']);
        saveas(fig_Te, filename_Te);
        fprintf("Teのグラフを %s に保存しました。\n", filename_Te);

        % neグラフのフルパスを生成して保存
        filename_ne = fullfile(saveDir, [base_filename, '_ne.png']);
        saveas(fig_ne, filename_ne);
        fprintf("neのグラフを %s に保存しました。\n", filename_ne);
    end
end