function [] = MP_plot(MP, MPdata2D, colorplot, date, savepath)
axisfontsize = 8;
switch colorplot
    case {'flow_t'}
        figure('Position', [0 0 1500 1500],'visible','on');
        for i = 1:MP.tate*MP.yoko
            offset_t = 1; % 何番目の要素からプロットを始めるか、基本的には1
            idx_t = offset_t + (i-1)*MP.dt*10;
            subplot(MP.tate,MP.yoko,i)
            switch colorplot
                case 'flow_t'
                    contourf(MPdata2D.zq,MPdata2D.rq,squeeze(MPdata2D.flow_forplot(idx_t,:,:)),100,'edgecolor','none');
                    c = colorbar;
                    clim([-40 40])
                    c.Label.String = 'Vt [km/s]';
                    colormap(redblue(300))
                    hold on
            end
            % フローを矢印で表示するなら必要そう
            % if MP.vector
            %     q = quiver(MPdata2D.zq,MPdata2D.rq,squeeze(MPdata2D.Ez(idx_t,:,:)),squeeze(MPdata2D.Er(idx_t,:,:)),3);
            %     q.Color = "k";
            % end
            % マッハプローブの位置をプロット
            [zprobe, rprobe] = meshgrid(MPdata2D.zprobe,MPdata2D.rprobe);
            scatter(zprobe, rprobe, 5, 'black', 'filled')
            hold on
            title([num2str(MPdata2D.trange(idx_t)+MP.delay) 'us'])
            xlabel('z [m]', 'FontSize', axisfontsize);
            ylabel('r [m]', 'FontSize', axisfontsize);
            xlim([-0.08 0.056])
            ylim([0.15 0.33])
            daspect([1 1 1])
        end
    otherwise
        return
end
min_shots_num = min(MP.shotlist);
max_shots_num = max(MP.shotlist);
mkdir(strcat(savepath, '\', num2str(date), '\output'));
saveas(gcf, strcat(savepath, '\', num2str(date), '\output\', num2str(min_shots_num), '-', num2str(max_shots_num), '_', num2str(MP.dt), '-us_mesh_', num2str(MP.mesh), '_mesh_r_', num2str(MP.meshr)), 'png');
end

