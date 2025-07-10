%%% cal flow with mean & sem %%%
function MPdata2D = MP_calc(pathname,MP, MPset)
n_ch = 9;%マッハプローブCH数

% 保存してあるものを使うならコメンならコメントアウト解除
% savename = [pathname,'\\mat\\mesh',num2str(MP.mesh),'_',num2str(MP.date),'_shot',num2str(MP.shotlist(1)),'-',num2str(MP.shotlist(end)),'-', num2str(MP.tate),'×',num2str(MP.yoko),'.mat'];
% if exist(savename,"file")
%     load(savename,'MPdata2D')
% else
    z = linspace(-0.08, 0.08, MP.mesh);%プロットメッシュZ座標[m]
    z_probe = [-0.08, -0.056, -0.032, -0.016, 0, 0.016, 0.032, 0.056, 0.08]; % マッハプローブ計測点Z座標[m]
    % ng_ch = [2, 3, 4, 5, 8, 9]; % 死んだCH 6/11
    % ng_ch = [3, 4, 5, 8, 9]; % 死んだCH 6/10
    ng_ch = [3, 4, 5, 7, 8, 9];
    z_probe(ng_ch) = [];
    MPdata2D.zprobe = z_probe;
    r = linspace(min(MP.rlist),max(MP.rlist),MP.meshr)*1E-3;%プロットメッシュR座標[m]
    r_probe = unique(MP.rlist)*1E-3; % マッハプローブ計測点R座標[m]
    MPdata2D.rprobe = r_probe;

    cnt_r = zeros(numel(r_probe),1);
    flow = zeros(numel(MP.trange),n_ch,numel(r_probe),MP.shot);
    for i = 1:numel(MP.shotlist)
        idx_r = find(r_probe==MP.rlist(i)*1E-3);
        cnt_r(idx_r) = cnt_r(idx_r) + 1;
        % filename = sprintf("%s%03d%s",[pathname.MP '\\' num2str(MP.date) '/ES_' num2str(MP.date)], MP.shotlist(i), '.csv');
        shotnum = MP.shotlist(i);
        if shotnum < 10
            shotnum = sprintf('%02d', shotnum);
        else 
            shotnum = num2str(shotnum);
        end
        filename = strcat(pathname, '\\',  num2str(MP.date), '\\I_ratio\\', shotnum, '.csv');
        MPdata = readmatrix(filename,'Range',sprintf('B%d:J%d',MP.trange(1)*10+2,MP.trange(end)*10+2));
        MPdata = MPdata ./ MPset.K * MPset.cs * 1E-3;
        flow(:,:,idx_r,cnt_r(idx_r)) = MPdata;
        
    end
    flow(:,ng_ch,:,:) = [];%死んだCHを除去
    MPdata2D.flow = flow;

    % 平均と標準誤差を格納するための変数を初期化
    flow_mean = zeros(size(MPdata2D.flow, 1), size(MPdata2D.flow, 2), size(MPdata2D.flow, 3));
    flow_error = zeros(size(MPdata2D.flow, 1), size(MPdata2D.flow, 2), size(MPdata2D.flow, 3));
    for id_R = 1:numel(r_probe)
        % 各rに対するshotの数を格納
        max_count = cnt_r(id_R);         
        % 平均を計算
        for count = 1:max_count
            data = MPdata2D.flow(:,:,id_R,count);
            flow_mean(:,:,id_R) = flow_mean(:,:,id_R) + data;
        end
        flow_mean(:,:,id_R) = flow_mean(:,:,id_R) / max_count;
        % 標準誤差を計算
        for count = 1:max_count
            data = MPdata2D.flow(:,:,id_R,count);
            flow_error(:,:,id_R) = flow_error(:,:,id_R) + (data - flow_mean(:,:,id_R)).^2;
        end
        flow_error(:,:,id_R) = sqrt(flow_error(:,:,id_R) / (max_count - 1))/sqrt(max_count);
    end
    MPdata2D.flow_mean = flow_mean;
    MPdata2D.flow_error = flow_error;

    [MPdata2D.zq,MPdata2D.rq] = meshgrid(z,r);
    MPdata2D.trange = MP.trange;
    MPdata2D.flow_forplot = zeros(numel(MP.trange),MP.meshr,MP.mesh);
    for i = 1:numel(MP.trange)
        MPdata2D.flow_forplot(i,:,:) = griddata(z_probe,r_probe,squeeze(flow_mean(i,:,:))',MPdata2D.zq,MPdata2D.rq);
    end
    savename = strcat(pathname, '\', num2str(MP.date),'\', shotnum, '_', num2str(MP.trange(1)*10), '-', num2str(MP.trange(end)*10), '.mat');
    save(savename,'MPdata2D')
% end
