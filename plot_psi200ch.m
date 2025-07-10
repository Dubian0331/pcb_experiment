function plot_psi200ch(PCB,pathname)
shot = PCB.shot;
trange = PCB.trange;
start = PCB.start;

[grid2D,data2D] = process_PCBdata_200ch(PCB,pathname);
% [grid2D,data2D] = process_PCBdata_200ch(PCB,pathname);

% ***********************************************

if isstruct(grid2D)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
    return
end

[magAxisList,xPointList] = get_axis_x_multi(grid2D,data2D); %時間ごとの磁気軸、X点を検索

% プロット部分
figure('Position', [0 0 1500 1500],'visible','on');

% figure('Position', [0 0 1500 1500],'visible','on');
% start=30;
dt = PCB.dt;
%  t_start=470+start;
Et_t = zeros(1,16);
E_eff = zeros(1,16);
mergingRatio = zeros(1,16);
% PCB.type = 0;
times = trange(1)+start:dt:trange(1)+start+dt*15;
 for m=1:16 %図示する時間
     i=start+(m-1)*dt; %start+400μsからdtごとに取得
     t=trange(i+1); %startが60なら460μsのデータを取る
     subplot(4,4,m);

     switch(PCB.type)
        case 0
            contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),10,'black','LineWidth',2);
        case 1
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none');clim([-5e-3,5e-3])%psi
        case 2
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none');clim([-0.1,0.1])%Bz
        case 3
            % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),40,'LineStyle','none');clim([0,0.3])%Bt
            % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt_th(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
        case 4
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Br(:,:,i),30,'LineStyle','none');clim([-0.1,0.1])%Br
        case 5
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none');clim([-500,400])%Et
            % 時刻iにおけるx点の[r,z]座標（インデックス）を取得
            idxR = knnsearch(grid2D.rq(:,1),xPointList.r(i));
            idxZ = knnsearch(grid2D.zq(1,:).',xPointList.z(i));
            % ±3程度の範囲でEtの平均を計算
            Et_t(m) = mean(data2D.Et(max(1,idxR-1):min(idxR+1,numel(grid2D.rq(:,1))),max(1,idxZ-1):min(numel(grid2D.zq(1,:)),idxZ+1),i),'all');
        case 6
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none');clim([-0.8*1e+6,0.8*1e+6]);%clim([-0.8*1e+6,0]) %jt%カラーバーの軸の範囲
        case 7
            Bp = sqrt(data2D.Bz(:,:,i).^2+data2D.Br(:,:,i).^2);
            % GFR = data2D.Bt_th(:,:,i)./Bp;
            % contourf(grid2D.zq(1,:),grid2D.rq(:,1),GFR,30,'LineStyle','none');clim([0,20]);%GFR
            % % 時刻iにおけるx点の[r,z]座標（インデックス）を取得
            % idxR = knnsearch(grid2D.rq(:,1),xPointList.r(i));
            % idxZ = knnsearch(grid2D.zq(1,:).',xPointList.z(i));
            % % ±3程度の範囲でGFRの平均（最大？）を計算
            % maxGFR = max(GFR(max(1,idxR-3):min(idxR+3,numel(grid2D.rq(:,1))),max(1,idxZ-3):min(numel(grid2D.zq(1,:)),idxZ+3)),[],'all');
            % disp(maxGFR);
            [zq,rq] = meshgrid(linspace(-0.1,0.1,200),linspace(0.2,0.32,200));
            Bp_q = griddata(grid2D.zq(1,:),grid2D.rq(:,1),Bp,zq,rq);
            Bt_th_q = griddata(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt_th(:,:,i),zq,rq);
            Et_q = griddata(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Et(:,:,i),zq,rq);
            GFR_q = abs(Bt_th_q./Bp_q);
            % contourf(zq(1,:),rq(:,1),GFR_q,50,'LineStyle','none');clim([0,50]);%GFR
            % disp(max(GFR_q,[],"all"))
            idxRq = knnsearch(rq(:,1),xPointList.r(i));
            idxZq = knnsearch(zq(1,:).',xPointList.z(i));
            E_eff_tmp = Et_q.*GFR_q;
            E_eff(m) = min(E_eff_tmp(max(1,idxRq-10):min(200,idxRq+10),max(1,idxZq-10):min(200,idxZq+10)),[],"all");
            contourf(zq(1,:),rq(:,1),-1*E_eff_tmp,50,'LineStyle','none');clim([1e3,Inf]);%GFR
            % Jt_q = griddata(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Jt(:,:,i),zq,rq);
            % eta_q = Et_q./Jt_q;
            % contourf(zq(1,:),rq(:,1),eta_q,linspace(-1e-2,1e-2,50),'LineStyle','none');clim([-1e-2,1e-2]);%resistivity
            % maxGFRq = max(GFR_q(max(1,idxRq-3):min(idxRq+3,numel(rq(:,1))),max(1,idxZq-3):min(numel(zq(1,:)),idxZq+3)),[],'all');
            % disp(maxGFRq);
        case 8
            Bp = sqrt(data2D.Bz(:,:,i).^2+data2D.Br(:,:,i).^2);
            % contourf(grid2D.zq(1,:),grid2D.rq(:,1),Bp,30,'LineStyle','none');clim([0 0.1]);%Bp
            [zq,rq] = meshgrid(linspace(-0.1,0.1,200),linspace(0.2,0.32,200));
            Bp_q = griddata(grid2D.zq(1,:),grid2D.rq(:,1),Bp,zq,rq);
            contourf(zq(1,:),rq(:,1),Bp_q,50,'LineStyle','none');clim([0,0.01]);%GFR
            disp(min(Bp_q,[],"all"));
        case 9
            B = sqrt(data2D.Bz(:,:,i).^2+data2D.Br(:,:,i).^2+data2D.Bt_th(:,:,i).^2);
            contourf(grid2D.zq(1,:),grid2D.rq(:,1),B,30,'LineStyle','none');clim([0 1]);%Bp  
     end
     mergingRatio(m) = xPointList.psi(i)/min(magAxisList.psi(:,i));

     % ポロイダル磁場ベースでのプロット
     % Bp = sqrt((data2D.Br(:,:,i)).^2+(data2D.Bz(:,:,i)).^2);
     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),Bp,30,'LineStyle','none');clim([0,0.02]);
     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),abs(data2D.Bt(:,:,i)./data2D.Bz(:,:,i)),30,'LineStyle','none');clim([-1,1]);
     % sepa = zeros(size(data2D.psi(:,:,i)));
     % x_psi = xPointList.psi(i);
     % sepa(((data2D.psi(:,:,i)<=x_psi+5e-4&data2D.psi(:,:,i)>=x_psi-5e-5)|(data2D.psi(:,:,i)<=x_psi*1.1&data2D.psi(:,:,i)>=x_psi*0.9)))=1;
     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),sepa,'LineStyle','none');
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none');clim([-5e-3,5e-3])%psi
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none');clim([-0.1,0.1])%Bz
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Br(:,:,i),30,'LineStyle','none');clim([-0.1,0.1])%Br
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),40,'LineStyle','none');clim([0,0.3])%Bt
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-100e-3:0.5e-3:100e-3,'LineStyle','none')
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none');clim([-0.8*1e+6,0]);%clim([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none');clim([-500,400])%Et
    % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Bt_th(:,:,i),20,'LineStyle','none');clim([0,0.3])%Bt
    colormap(jet)
    axis image
    axis tight manual
%     caxis([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
    % clim([-0.1,0.1])%Bz
     % clim([0,0.3])%Bt
    % clim([-5e-3,5e-3])%psi
    % clim([-500,400])%Et
%     colorbar('Location','eastoutside')
    %カラーバーのラベル付け
%     c = colorbar;
%     c.Label.String = 'Jt [A/m^{2}]';
    hold on
%     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
%     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')

    contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')

    % contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.2e-3:40e-3],'black','LineWidth',1)
%     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
%     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
     % plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置

    % timeIndex = find(trange==t);
    % [magaxis,xpoint] = get_axis_x(grid2D,data2D,t);
    % plot(magaxis.z,magaxis.r,'ro');
    % plot(xpoint.z,xpoint.r,'rx');

    plot(magAxisList.z(:,i),magAxisList.r(:,i),'ko');
    plot(xPointList.z(i),xPointList.r(i),'kx');

    if PCB.type == 7 || PCB.type == 8
        xlim([-0.1,0.1]);ylim([0.2,0.32]);
    end
    hold off
    title(string(t)+' us')
%     xlabel('z [m]')
%     ylabel('r [m]')
 end

 sgtitle(strcat('shot',num2str(shot)));

 if PCB.type==5
     figure;plot(times,Et_t,'LineWidth',3);
     ax = gca;
     ax.FontSize = 18;
     xlabel('time [us]');ylabel('Reconnection electric field [V/m]');
 elseif PCB.type == 7
     figure;plot(times,E_eff,'LineWidth',3);
    xlabel('time [us]');
    %  figure;plot(mergingRatio,E_eff,'LineWidth',3);
    %  ax = gca;
    %  ax.FontSize = 18;
    % %  xlim([460 468]);
    %  xlabel('Merging ratio');
 end
 % drawnow

end










% %%%%%%%%%%%%%%%%%%%%%%%%
% %以下、local関数
% %%%%%%%%%%%%%%%%%%%%%%%%

% function plot_psi200ch(date, shot, tfshot, pathname, n_r ,n_z, i_EF,trange,doCheck)

% %計算有無
% % def) filename: processeddataのファイル名
% filename = strcat(pathname.processeddata,'\a039_',num2str(shot(1)),'.mat');

% %較正係数のバージョンを日付で判別
% sheets = sheetnames('coeff280ch.xlsx');
% sheets = str2double(sheets);
% sheet_date=max(sheets(sheets<=date));
% C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
% r_shift_pcb = 0.000;
% ok = logical(C(:,14));
% dtacq_num_list = C(:,1);
% dtaq_ch = C(:,2);
% polarity=C(:,13);
% coeff=C(:,12);
% zpos=C(:,9);
% rpos=C(:,10)+r_shift_pcb;
% ch=C(:,7);

% if ismember(39,dtacq_num_list)
%     filename1 = strcat(pathname.rawdata,'\rawdata_dtacq',num2str(39),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
%     if exist(filename1,"file")==0
%         disp(['File:',filename1,' does not exit']);
%         disp('No rawdata file of a039 -- Start generating!')
%         rawdataPath = pathname.rawdata;
%         save_dtacq_data(39, shot(1), tfshot(1),rawdataPath)
%         % doCalculation=true;
%         % return
%     end
%     a039_raw = importdata(filename1);
% end
% if ismember(40,dtacq_num_list)
%     filename2 = strcat(pathname.rawdata,'/rawdata_dtacq',num2str(40),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
%     if exist(filename2,"file")==0
%         disp(['File:',filename2,' does not exit']);
%         % return
%         % doCalculation=true;
%         disp('No rawdata file of a040 -- Start generating!')
%         rawdataPath = pathname.rawdata;
%         save_dtacq_data(40, shot(2), tfshot(2),rawdataPath)
%     end
%     a040_raw = importdata(filename2);
% end

% raw = zeros(1000,length(dtaq_ch));
% for i = 1:length(dtaq_ch)
%     if dtacq_num_list(i) == 39
%         raw(:,i) = a039_raw(:,dtaq_ch(i));
%     elseif dtacq_num_list(i) == 40
%         raw(:,i) = a040_raw(:,dtaq_ch(i));
%     end
% end

% b=raw.*coeff';%較正係数RC/NS
% b=b.*polarity';%極性揃え

% %デジタイザchからプローブ通し番号順への変換
% bz=zeros(1000,100);
% bt=bz;
% ok_bz=false(100,1);
% ok_bt=ok_bz;
% zpos_bz=zeros(100,1);
% rpos_bz=zpos_bz;
% zpos_bt=zpos_bz;
% rpos_bt=zpos_bz;

% %digital filter
% windowSize = 3;
% bb = (1/windowSize)*ones(1,windowSize);
% aa = 1;

% for i=1:length(ch)
%     b(:,i) = filter(bb,aa,b(:,i));
%     b(:,i) = b(:,i) - mean(b(1:40,i));
%     if rem(ch(i),2)==1
%         bz(:,ceil(ch(i)/2))=b(:,i);
%         ok_bz(ceil(ch(i)/2))=ok(i);
%         zpos_bz(ceil(ch(i)/2))=zpos(i);
%         rpos_bz(ceil(ch(i)/2))=rpos(i);
%     elseif rem(ch(i),2)==0
%         bt(:,ch(i)/2)=b(:,i);
%         ok_bt(ceil(ch(i)/2))=ok(i);
%         zpos_bt(ceil(ch(i)/2))=zpos(i);
%         rpos_bt(ceil(ch(i)/2))=rpos(i);
%     end
% end

% % zprobepcb    = [-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17];
% zprobepcb    = [-0.2975,-0.255,-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17,0.255,0.2975];
% rprobepcb    = [0.06,0.09,0.12,0.15,0.18,0.21,0.24,0.27,0.30,0.33]+r_shift_pcb;
% rprobepcb_t  = [0.07,0.10,0.13,0.16,0.19,0.22,0.25,0.28,0.31,0.34]+r_shift_pcb;
% % rprobepcb_new = [0.18,0.195,0.21,0.225,0.24,0.255,0.27,0.285,0.3,0.315];
% rprobepcb_new = [0.18,0.195,0.21,0.225,0.24,0.255,0.27,0.285,0.3,0.315]-0.12;
% rprobepcb_t_new = rprobepcb_new+0.01;
% [zq_probepcb_bt,rq_probepcb_bt]=meshgrid(zprobepcb,rprobepcb_t);
% [zq_probepcb,rq_probepcb]=meshgrid(zprobepcb,rprobepcb);
% ok_bt_matrix = false(length(rprobepcb),length(zprobepcb));
% ok_bz_matrix = false(length(rprobepcb),length(zprobepcb));
% for i = 1:length(ok_bt)
%     if rpos_bt(i) > (r_shift_pcb)
%         index_r = (abs(rpos_bt(i)-rprobepcb_t)<0.001);index_z = (zpos_bt(i)==zprobepcb);
%         ok_bt_matrix = ok_bt_matrix + rot90(index_r,-1)*index_z*ok_bt(i);
%     end
%     index_r = (abs(rpos_bz(i)-rprobepcb)<0.001);index_z = (zpos_bz(i)==zprobepcb);
%     ok_bz_matrix = ok_bz_matrix + rot90(index_r,-1)*index_z*ok_bz(i);
% end

% grid2D_pcb=struct(...
%     'zq',zq_probepcb,...
%     'rq',rq_probepcb,...
%     'zprobepcb',zprobepcb,...
%     'rprobepcb',rprobepcb,...
%     'rprobepcb_t',rprobepcb_t,...
%     'rprobepcb_new',rprobepcb_new,...
%     'rprobepcb_t_new',rprobepcb_t_new,...
%     'ok_bz_matrix',ok_bz_matrix,...
%     'ok_bt_matrix',ok_bt_matrix,...
%     'zq_probepcb_bt',zq_probepcb_bt,...
%     'rq_probepcb_bt',rq_probepcb_bt);
% clear zq_probepcb rq_probepcb ok_bz_matrix ok_bt_matrix zq_probepcb_bt rq_probepcb_bt


% % ******************* old way ********************
% for i=1:size(trange,2)
%     t=trange(i);

%     %Bzの二次元補間(線形fit)
%     B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D_pcb, bt, ok_bt, t);

% end
% % ***********************************************

% B_z_combined = bz;
% B_t_combined = bt;
% [grid2D_new_zq,grid2D_new_rq] = meshgrid(linspace(min(grid2D_pcb.zq(1,:)),max(grid2D_pcb.zq(1,:)),n_z),...
%     linspace(min(grid2D_pcb.rq(:,1)),max(grid2D_pcb.rq(:,1)),n_r));
% grid2D = struct(...
%     'zq',grid2D_new_zq,...
%     'rq',grid2D_new_rq,...
%     'zprobepcb',grid2D_pcb.zprobepcb,...
%     'rprobepcb',grid2D_pcb.rprobepcb,...
%     'zq_pcb_bt',grid2D_pcb.zq_probepcb_bt,...
%     'rq_pcb_bt',grid2D_pcb.rq_probepcb_bt,...
%     'okbz_pcb',grid2D_pcb.ok_bz_matrix, ...
%     'okbt_pcb',grid2D_pcb.ok_bt_matrix);
% data2D =struct(...
%     'psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%     'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%     'Bt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%     'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%     'Jz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%     'Jr',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%     'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%     'P',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%     'Lambda',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),...
%     'trange',trange);

% %プローブチェック
% if doCheck
%     probecheck_script(ok_bz,ok_bt,bz,bt);
% end

% %data2Dcalc.m
% r_EF   = 0.5 ;
% n_EF   = 234. ;
% if date<221119
%     z1_EF   = 0.68;
%     z2_EF   = -0.68;
% else
%     z1_EF   = 0.78;
%     z2_EF   = -0.78;
% end
% [Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
% clear EF r_EF n_EF i_EF z_EF

% for i=1:size(trange,2)
%     t=trange(i);
%     %%Bzの二次元補間(線形fit)
%     vq = bz_rbfinterp(rpos_bz, zpos_bz, grid2D, B_z_combined, ok_bz, t);
%     B_z = -Bz_EF+vq;
%     B_t = bz_rbfinterp(rpos_bt, zpos_bt, grid2D, B_t_combined, ok_bt, t);
    
%     %PSI計算
%     data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),2*pi*B_z.*grid2D.rq(:,1),1);
    
%     %このままだと1/2πrが計算されてないので
%     [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
%     data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
% %     data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
%     data2D.Bz(:,:,i)=B_z;
%     data2D.Bt(:,:,i)=B_t;
%     [curlt,~] = curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i));
%     data2D.Jt(:,:,i)= curlt/(4*pi*1e-7);
%     [~,dRBt_dR] = gradient(grid2D.rq.*data2D.Bt(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1));
%     [dBt_dZ,~] = gradient(data2D.Bt(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1));
%     data2D.Jz(:,:,i) = 1./grid2D.rq.*dRBt_dR./(4*pi*1e-7);
%     data2D.Jr(:,:,i) = -dBt_dZ./(4*pi*1e-7);
%     data2D.Lambda(:,:,i) = (2*pi*grid2D.rq.*data2D.Bt(:,:,i))./(data2D.psi(:,:,i));
%     data2D.P(:,:,i) = cumtrapz(grid2D.rq(:,1),data2D.Jt(:,:,i).*data2D.Bz(:,:,i)-data2D.Jz(:,:,i).*data2D.Bt(:,:,i),1);
% end

% if isstruct(grid2D_pcb)==0 %もしdtacqデータがない場合次のloopへ(データがない場合NaNを返しているため)
%     return
% end



% %------------------------------------%
% % マッハプローブ用 %
% mpz = [-0.08, -0.056, -0.032, -0.016, 0, 0.016, 0.032, 0.056, 0.08];
% mpz = [-0.08, -0.056, -0.032, 0, -0.016, 0, 0.032, 0.056];
% mpr = 0.15:0.03:0.33;

% % 90度回転ver
% mpz = 0.08;
% probe_length = [0.08, 0.056, 0.032, 0.016, 0, -0.016, -0.032, -0.056, -0.08];
% probe_length = probe_length + 0.08;
% position = 0.4; % ここを挿す位置によって調整
% position = 0.525 + 0.02 - position;
% mpr = sqrt(ones(1, 9) * position.^2 + probe_length.^2);
% %------------------------------------%



% figure('Position', [0 0 1500 1500],'visible','on');
% start=60;
% dt = 3;
% %  t_start=470+start;
%  for m=1:16 %図示する時間
%      i=start+m.*dt; %end
%      t=trange(i);
%      subplot(4,4,m)
%     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),30,'LineStyle','none')
%     contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.psi(:,:,i),40,'LineStyle','none')
%     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bt(:,:,i),-50e-3:0.2e-3:50e-3,'LineStyle','none')
%     % contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Jt(:,:,i),30,'LineStyle','none')
% %     contourf(grid2D.zq(1,:),grid2D.rq(:,1),-1.*data2D.Et(:,:,i),20,'LineStyle','none')
%     colormap(jet)
%     axis image
%     axis tight manual
% %     caxis([-0.8*1e+6,0.8*1e+6]) %jt%カラーバーの軸の範囲
%     % caxis([-0.05,0.05])%Bz
%     % clim([-0.05,0.05])%Bt
%     clim([-8e-3,8e-3]);%psi
% %     caxis([-500,400])%Et
%     % colorbar('Location','eastoutside')
%     %カラーバーのラベル付け
%     c = colorbar;
%     % c.Label.String = 'Bt [T]';
%     c.Label.FontSize = 12;
%     hold on
% %     plot(grid2D.zq(1,squeeze(mid(:,:,i))),grid2D.rq(:,1))
% %     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
% %     contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),20,'black')
%      contour(grid2D.zq(1,:),grid2D.rq(:,1),squeeze(data2D.psi(:,:,i)),[-20e-3:0.4e-3:40e-3],'black','LineWidth',1)
% %     plot(grid2D.zq(1,squeeze(mid(opoint(:,:,i),:,i))),grid2D.rq(opoint(:,:,i),1),"bo")
% %     plot(grid2D.zq(1,squeeze(mid(xpoint(:,:,i),:,i))),grid2D.rq(xpoint(:,:,i),1),"bx")
%      % plot(ok_z,ok_r,"k.",'MarkerSize', 6)%測定位置
%     % for r = mpr
%     %     scatter(mpz, r, 10, 'black', 'filled') %マッハプローブ用
%     % end
%     hold off
%     title(string(t)+' us')
%     xlim([-0.17 0.17]);
%     % ylim([0.06, 0.34]);
%     xlabel('z [m]', 'FontSize',12)
%     ylabel('r [m]','FontSize',12)
%  end
% % figure('Position', [0 0 1500 1500],'visible','on');
% % for m=1:16 %図示する時間
% %     i=start+m.*dt; %end
% %     t=trange(i);
% %     subplot(4,4,m)
% %     plot(grid2D.zq(1,:),data2D.Bt(25,:,i))
% %     title(string(t)+' us')
% %     ylim([-0.03 0.03])
% % end
% % saveas(gcf,strcat('/Users/yunhancai/Downloads/files/a039_',num2str(shot)))
% % saveas(gcf,strcat('/Users/yunhancai/Downloads/a039_',num2str(shot),'.png'))
% % close

% end