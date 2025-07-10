function vq = b_interp(rpos, zpos, grid2D, bz, ok, t)
%磁場データの線形補間
%   input:bz[samplingnumber×ch],rpos(ch),zpos(ch),grid2D{rq(n×n),zq(n×n)},ok(128(true or false)),t[us]
%   output:vq(n×n) ：bzがn×nのグリッドにrbf補間されたもの

%%無視するチャンネルを除いたbzの散布データ（okのチャンネルのみ残す）
x=zpos(ok);
y=rpos(ok);
z=double(bz(t,ok))';

f = fit([x,y],z,'lowess','Span',0.1);%  f = fit([x,y],z,'cubicinterp');
vq=f(grid2D.zq,grid2D.rq);
% plot(f,[x,y],z,'Style','Residuals'); %vq=f(grid2D.zq,grid2D.rq);;contourf(vq)
vq = fillmissing(vq,'linear');

% %%rbf補間
% %%スムージングと関数の選択
% smoothval=0.05;
% func='multiquadric';%Gaussian, Linear, Cubic,multiquadric, Thinplate から選べる
% const=0.015;
% vq= rbfinterp([grid2D.zq(:)'; grid2D.rq(:)'], rbfcreate([x' ; y'], z','RBFFunction',func,'RBFConstant',const,'RBFSmooth', smoothval ));%,'RBFSmooth', smoothval,'RBFConstant',const));
% vq = reshape(vq, size(grid2D.zq));
end