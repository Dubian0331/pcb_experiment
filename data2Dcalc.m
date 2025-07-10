function data2D = data2Dcalc(EF, grid2D, n, trange, rpos, zpos, bz, ok)
r_EF   = 0.5 ;
n_EF   = 234. ;
i_EF    =EF;
z1_EF   = 0.68;%0.702
z2_EF   = -0.68;
[Bz_EF,~] = B_EF(z1_EF,z2_EF,r_EF,i_EF,n_EF,grid2D.rq,grid2D.zq,false);
clear EF r_EF n_EF i_EF z_EF
%% 
% psi計算用 (Bzの二次元補間→PSI計算)

%trange=(468:1:477);
data2D=struct('psi',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Bz',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Br',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Jt',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'Et',zeros(size(grid2D.rq,1),size(grid2D.rq,2),size(trange,2)),'trange',trange);
minpsi=zeros(size(grid2D.rq,1),size(trange,2));
minind=zeros(size(grid2D.rq,1),size(trange,2));
for i=1:size(trange,2)
    t=trange(i);
    %%Bzの二次元補間（rbfinterp）
    vq = bz_rbfinterp(rpos, zpos, grid2D, bz, ok, t);
    
    B_z = -Bz_EF+vq;
    %%PSI計算
    data2D.psi(:,:,i) = cumtrapz(grid2D.rq(:,1),B_z,1);
    %このままだと1/2πrが計算されてないので
    [data2D.Br(:,:,i),data2D.Bz(:,:,i)]=gradient(data2D.psi(:,:,i),grid2D.zq(1,:),grid2D.rq(:,1)) ;
    data2D.Br(:,:,i)=-data2D.Br(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Bz(:,:,i)=data2D.Bz(:,:,i)./(2.*pi.*grid2D.rq);
    data2D.Jt(:,:,i)= curl(grid2D.zq(1,:),grid2D.rq(:,1),data2D.Bz(:,:,i),data2D.Br(:,:,i))./(4*pi*1e-7);
end
data2D.Et=diff(data2D.psi,1,3).*1e+6; 
%diffは単なる差分なので時間方向のsizeが1小さくなる %ステップサイズは1us
data2D.Et=data2D.Et./(2.*pi.*grid2D.rq);
%data2D=struct('psi',psi,'Bz',Bz,'Br',Br,'Jt',Jt,'trange',trange);

clear vq bz Bz_EF B_z ok ng rpos zpos i t
end