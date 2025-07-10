%%%%%%%%%%%pcb計測の二次元データdata2D（Bz,Br，Et，Jt）とグリッドgrid2D（zq,rq）が出てくる。
%%%% input （日付、データックのショット番号、差し引きショット、見たい時間、nはグリッドの個数、EFの電流）
%%%%%%↑　これはtestの中にある。psisave.mlx内のfunction [date, shot, TF_shot, offset_TF,i_EF,start,Doppler_t] = getinput(T,IDX)を関数化して呼び出すのが楽
function [grid2D, data2D] = pcbdata(date, d_tacq,d_tacqTF,trange, coeff, n,EF)
% % mdsplusを通じて読み込み
% [rawdata]=getvalue(d_tacq,d_tacqTF);

%localに保管したrawdataから読み込み
pathname.rawdata='C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\rawdata_a038\'; %rawdataの保管場所
filename=strcat(pathname.rawdata,'rawdata_dtacq',num2str(d_tacq),'.mat');

% %noTF
% pathname.rawdata='C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\rawdata_a038_noTF\'; %rawdataの保管場所
% filename=strcat(pathname.rawdata,'rawdata_noTF_dtacq',num2str(d_tacq),'.mat');

load(filename,'rawdata');

if numel(rawdata)< 500
    grid2D=NaN;
    data2D=NaN;
    return
end
if exist(coeff)==0
    load('rc_coeff2020.mat')
end     
[ok, bz, rpos, zpos,p_ch] = getpcbbz(rawdata, coeff,date);

% チャンネルごとの生信号のプロット
bz=smoothdata(bz,1);%各ch内のデータを時間に関してsmoothing

r = 5;%プローブ本数＝グラフ出力時の縦に並べる個数
col1 = 12;%1枚目のグラフ出力時の横に並べる個数
col2 = 13;%2枚目のグラフ出力時の横に並べる個数
y_upper_lim = 0.03;%2.5e-3;%0.1;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.03;%2e-3;%-0.1;%縦軸プロット領域（b_z下限）
t_start=330;%455;%横軸プロット領域（開始時間）
t_end=600;%横軸プロット領域（終了時間）
r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数
plotbzsignal(y_upper_lim, col2, col1, t_end, p_ch, y_lower_lim, t_start, bz, ok, r_ch, r);
clear r col1 col2 y_upper_lim y_lower_lim t t_start t_end

[zq,rq]=meshgrid(linspace(min(zpos),max(zpos),n),linspace(min(rpos),max(rpos),n));
grid2D=struct('zq',zq,'rq',rq);
clear zq rq

data2D = data2Dcalc(EF, grid2D, n, trange, rpos, zpos, bz, ok);
end
