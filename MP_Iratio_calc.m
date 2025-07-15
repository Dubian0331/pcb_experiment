function MPdata2D = MP_Iratio_calc(pathname,MP, MPset)
    n_ch = 9;
    z = linspace(-0.08, 0.08, MP.mesh);
    z_probe = [-0.08, -0.056, -0.032, -0.016, 0, 0.016, 0.032, 0.056, 0.08];
    ng_ch = MP.ng_ch;
    z_probe(ng_ch) = [];
    MPdata2D.zprobe = z_probe;
    r = linspace(min(MP.rlist),max(MP.rlist),MP.meshr)*1E-3;
    r_probe = unique(MP.rlist)*1E-3;
    MPdata2D.rprobe = r_probe;
    cnt_r = zeros(numel(r_probe),1);

    I_ratio_4D = zeros(numel(MP.trange),n_ch,numel(r_probe),MP.shot);
    
    modified_folder_path = fullfile(pathname, num2str(MP.date), 'I_ratio_modified');
    if isfolder(modified_folder_path)
        target_folder = 'I_ratio_modified';
    else
        target_folder = 'I_ratio';
    end
    disp(['Reading from folder: ', target_folder]);

    for i = 1:numel(MP.shotlist)
        idx_r = find(r_probe==MP.rlist(i)*1E-3);
        cnt_r(idx_r) = cnt_r(idx_r) + 1;
        shotnum = MP.shotlist(i);
        if shotnum < 10
            shotnum = sprintf('%02d', shotnum);
        else 
            shotnum = num2str(shotnum);
        end
        filename = strcat(pathname, '\\',  num2str(MP.date), '\\', target_folder, '\\', shotnum, '.csv');
        MPdata = readmatrix(filename,'Range',sprintf('B%d:J%d',MP.trange(1)*10+2,MP.trange(end)*10+2));
        
        % Flow速度への変換式をコメントアウト
        % MPdata = MPdata ./ MPset.K * MPset.cs * 1E-3;
        
        I_ratio_4D(:,:,idx_r,cnt_r(idx_r)) = MPdata;
    end

    I_ratio_4D(:,ng_ch,:,:) = []; %死んだCHを除去
    
    I_ratio_mean = zeros(size(I_ratio_4D, 1), size(I_ratio_4D, 2), size(I_ratio_4D, 3));
    I_ratio_error = zeros(size(I_ratio_4D, 1), size(I_ratio_4D, 2), size(I_ratio_4D, 3));

    for id_R = 1:numel(r_probe)
        max_count = cnt_r(id_R);         
        % 平均を計算
        for count = 1:max_count
            data = I_ratio_4D(:,:,id_R,count);
            I_ratio_mean(:,:,id_R) = I_ratio_mean(:,:,id_R) + data;
        end
        I_ratio_mean(:,:,id_R) = I_ratio_mean(:,:,id_R) / max_count;
        % 標準誤差を計算
        for count = 1:max_count
            data = I_ratio_4D(:,:,id_R,count);
            I_ratio_error(:,:,id_R) = I_ratio_error(:,:,id_R) + (data - I_ratio_mean(:,:,id_R)).^2;
        end
        I_ratio_error(:,:,id_R) = sqrt(I_ratio_error(:,:,id_R) / (max_count - 1))/sqrt(max_count);
    end

    MPdata2D.I_ratio_raw = I_ratio_4D;
    MPdata2D.I_ratio_mean = I_ratio_mean;
    MPdata2D.I_ratio_error = I_ratio_error;
    
    [MPdata2D.zq,MPdata2D.rq] = meshgrid(z,r);
    MPdata2D.trange = MP.trange;

    I_ratio_forplot = zeros(numel(MP.trange),MP.meshr,MP.mesh);
    for i = 1:numel(MP.trange)
        I_ratio_forplot(i,:,:) = griddata(z_probe,r_probe,squeeze(I_ratio_mean(i,:,:))',MPdata2D.zq,MPdata2D.rq);
    end
    MPdata2D.I_ratio_forplot = I_ratio_forplot; % グリッドデータも名前を変更
    
    % ---matファイルを保存---
    save_folder = fullfile(pathname, 'mat', num2str(MP.date));
    if ~exist(save_folder, 'dir')
        mkdir(save_folder);
    end
    savename = fullfile(save_folder, [num2str(MP.shotlist(1)), '-', num2str(MP.shotlist(end)), '_Iratio.mat']);
    save(savename, 'MPdata2D');
end