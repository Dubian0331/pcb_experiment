function [] = MP_plot_with_psi(MP, MPdata2D, colorplot, date, data2D, grid2D, savepath)
    axisfontsize = 14;
    clim_values = [-40 40];  % カラーマップの範囲を指定
    
    % タイルレイアウトを作成
    % t = tiledlayout(MP.tate, MP.yoko, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    switch colorplot
        case {'flow_t'}
            figure('Position', [0 0 800 800],'visible','on', WindowState='maximized');
            tiledlayout(MP.tate, MP.yoko, 'TileSpacing', 'compact', 'Padding', 'compact');
            for i = 1:MP.tate * MP.yoko
                [~, offset_t] = min(abs(MP.trange - MP.plotstart)); % 指定した値からプロットを始める
                offset_mag_t = 399;
                idx_t = offset_t + (i - 1) * MP.dt * 10;
                nexttile
                hold on
                contourf(MPdata2D.zq, MPdata2D.rq, squeeze(MPdata2D.flow_forplot(idx_t, :, :)), 100, 'edgecolor', 'none');
                contourf(grid2D.zq, grid2D.rq, data2D.psi(:, :, MP.trange(idx_t) - offset_mag_t), -6e-3:0.4e-3:6e-3, 'k', 'Fill', 'off', 'LineWidth', 1);
                colormap(redblue(300))
                clim(clim_values)
                disp([num2str(data2D.trange(MP.trange(idx_t) - 399)), ', ', num2str(MPdata2D.trange(idx_t))]);
                
                % マッハプローブの位置をプロット
                % zprobe = [-0.08, -0.056, -0.032, -0.016, 0, 0.016, 0.032, 0.056, 0.08];
                % [zprobe, rprobe] = meshgrid(zprobe, MPdata2D.rprobe);
                % scatter(zprobe, rprobe, 20, 'black', 'filled')
                hold on
                title(['t=', num2str(MPdata2D.trange(idx_t) + MP.delay),' us'], 'FontSize', 14)

                % X軸とY軸のラベルを設定
                if mod(i, MP.yoko) == 1
                    ylabel('r [m]', 'FontSize', axisfontsize);
                    ax = gca;
                    ax.YAxis.FontSize = axisfontsize;
                else
                    set(gca, 'YTickLabel', []);
                end
                if i > MP.yoko * (MP.tate - 1)
                    ax = gca;
                    ax.XAxis.FontSize = axisfontsize;
                    xticks([-0.08 -0.032 0 0.032 0.08])
                    xlabel('z [m]', 'FontSize', axisfontsize);
                else
                    set(gca, 'XTickLabel', []);
                end

                xlim([-0.08 0.08])
                ylim([0.10 0.30])
                daspect([1 1 1])
            end
            % 共通のカラーバーを設定
            cb = colorbar;
            cb.Layout.Tile = 'east';
            cb.FontSize = 14;
            cb.Label.String = 'Vt [km/s]';
            drawnow;
        otherwise
            return
    end
    min_shots_num = min(MP.shotlist);
    max_shots_num = max(MP.shotlist);

    % --- 保存 ---
    mkdir(strcat(savepath, '\figure\', num2str(date)));
    saveas(gcf, strcat(savepath, '\figure\', num2str(date), '\', num2str(min_shots_num), '-', num2str(max_shots_num), '_', num2str(MP.dt), 'us_', num2str(MP.plotstart) ,'us-_mesh_', num2str(MP.mesh), '_mesh_r_', num2str(MP.meshr)), 'png');
end
