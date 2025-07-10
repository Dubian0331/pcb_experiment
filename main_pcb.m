clear all

for i = 1:36
%     load(['/Users/yunhancai/Google Drive/Data/pcb/pre_processed/a039_',num2str(i),'.mat']);
%     plot_psi_multi_pcb(data2D,grid2D,false,true,i);
%     merging_rate2(data2D.psi,grid2D.zq(1,:),rot90(grid2D.rq(:,1),3),400,600);
%     helicity_calculation(data2D,rot90(grid2D.rq(:,1),3),grid2D.zq(1,:),rot90(grid2D.rq(:,1),3),grid2D.zq(1,:));
%     lambda_calculation(data2D,grid2D)2
    rogowski(230412,10,0,i);
    saveas(gcf,strcat('/Users/yunhancai/Downloads/rogo_230412_',num2str(i),'.png'))
%     saveas(gcf,strcat('/Users/yunhancai/Google Drive/Data/pcb/a039_',num2str(i),'.png'))
%     disp(i)
    clearvars -except i
end

%plot_B_z_in_time(B_z,ch_dist,350,600)
%plot_B_z_vs_z(B_z_new,440,r_bz_probe,z_bz_probe,ch_dist);
%[psi_1,psi_2,psi_common] = merging_rate2(data2D.psi,grid2D.zq(1,:),rot90(grid2D.rq(:,1),3),450-300,600-300);
%plot_psi_at_t(Bzcoil,r_bz_probe,z_bz_probe,1,false,true,true,false)
%plot_B_z_in_rz(B_z_new,r_bz_probe,z_bz_probe);
%plot_B_z_in_rz(B_t_new,r_bt_probe,z_bt_probe);