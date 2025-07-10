function [] = movie_MP(pathname, MP,MPdata2D,colorplot)
sizen= 512; % よく分からないけど、このサイズが小さいと権限がないと言われるらしい
axisfontsize = 20;
titlefontsize = 20;
delaytime = 0.5; % 画像を送る間隔 (s)
switch colorplot
    case {'flow_t'}
        h = figure('Position', [0 0 1500 1500],'visible','on');
        % h.WindowState = 'maximized';
        % h = figure(WindowState="maximized");
        axis tight manual
        filename = strcat(pathname, '\\toroidal_flow.gif');
        for i = 1:30
            offset_t = 1; % 何番目の要素からプロットを始めるか、基本的には1
            dt = 1; % 1秒ごとにプロット
            idx_t = offset_t + (i-1)*dt*10;
            switch colorplot
                case 'flow_t'
                    contourf(MPdata2D.zq,MPdata2D.rq,squeeze(MPdata2D.flow_forplot(idx_t,:,:)),100,'edgecolor','none');
                    c = colorbar;
                    clim([-20 20])
                    c.Label.String = 'Vt [km/s]';
                    c.Label.FontSize = 20;
                    colormap(redblue(300))
                    hold on
            end
        hold on
        [zprobe, rprobe] = meshgrid(MPdata2D.zprobe,MPdata2D.rprobe);
        scatter(zprobe, rprobe, 5, 'black', 'filled')
        hold on
        title([num2str(MPdata2D.trange(idx_t)+MP.delay) 'us'], FontSize=titlefontsize)
        xlabel('z [m]', 'FontSize', axisfontsize);
        ylabel('r [m]', 'FontSize', axisfontsize);
        xlim([-0.08 0.056])
        ylim([0.15 0.33])
        drawnow
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,sizen);
        % Write to the GIF File
        if i == 1
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',delaytime);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',delaytime);
        end
        end
    otherwise
        return
end

end

