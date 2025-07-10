function [] = plot_MP_r_flow(MP, MPdata2D, ch)
axisfontsize = 8;
titlefontsize = 10;
titlefontsize2 = 12;

z_probe_original = [-0.08, -0.056, -0.032, -0.016, 0, 0.016, 0.032, 0.056, 0.08];
z = z_probe_original(ch);
id_z = find(MPdata2D.zprobe == z);
graphnum = 4; % 表示するグラフの数
start_id_t = 15; % MP.trange(1) + start_id_t + MP.delay [s] から描画
figure(WindowState="maximized")
for i = 1:graphnum
    subplot(1,graphnum,i)
    id_t = 1 + start_id_t*10 + (i-1)*MP.dt2*10;
    errorbar(MPdata2D.rprobe, squeeze(MPdata2D.flow_mean(id_t,id_z,:)), squeeze(MPdata2D.flow_error(id_t,id_z,:)))
    title([num2str(MPdata2D.trange(id_t)+MP.delay) 'us'], FontSize=titlefontsize);
    xlim([MPdata2D.rprobe(1) - 0.02, MPdata2D.rprobe(end) + 0.02]);
    ylim([-20, 20]);
    xlabel('r [m]', 'FontSize', axisfontsize);
    ylabel('Vt [km/s]', 'FontSize', axisfontsize);
    yline(0, 'Color', 'k', 'LineStyle', '--');% y = 0
    grid on;
end
sgtitle(strcat("z = ", num2str(z), " m"), fontsize=titlefontsize2);
