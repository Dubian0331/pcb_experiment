%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% マッハプローブ
% 流速をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% addpath '/Users/rsomeya/Documents/lab/matlab/common';
% run define_path.m
clc

% データとして採用するショット番号: case-I, 240610
MP.date = 240610;%【input】マッハプローブ計測日
MP.shotlist = [18, 20, 23, 26:27, 31:33, 35:38, 40:41, 44, 47];

% データとして採用するショット番号: case-O, 240611
% MP.date = 240611;%【input】マッハプローブ計測日
% MP.shotlist = [47, 48, 51:52, 54:60, 62, 65:70];

MP.delay = 0;  % トリガタイミングの時間ずれ（オシロスコープ）
MP.start = 480 + MP.delay;  %【input】プロット開始時刻[us]
MP.shot = 3;  % 同じ条件で打ったショットの最大値
MP.scandt = 0.1;  % 計算の刻み幅

% 測定データが入っているフォルダ
machpathname = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\\MachProbe_data';


MP.mesh = 100;  %【input】マッハプローブ補間メッシュ数（z方向）
MP.meshr = 100;  %【input】マッハプローブ補間メッシュ数（r方向）
MP.start = 470;
MP.end = 521;
MP.trange = MP.start:MP.scandt:MP.end;  %【input】計算時間範囲
MP.tate = 3;  %【input】プロット枚数(縦)
MP.yoko = 6;  %【input】プロット枚数(横)
MP.dt = 3;  %【input】プロット時間間隔[us] plot_MP用
MP.dt2 = 5;  %【input】プロット時間間隔[us] plot_MP_r_flow用

% MP.ng_ch = [3, 4, 5, 7, 8, 9];
% MP.ng_ch = [2, 3, 4, 5, 8, 9]; % 死んだCH 6/11
% MP.ng_ch = [3, 4, 5, 8, 9]; % 死んだCH 6/10
MP.ng_ch = [4, 5, 7, 8, 9];
MP.Vd2 = 40;  %【input】langmiur用の値（手前と書かれた方のバイアス電圧値）
MP.Vd3 = 40;  %【input】langmiur用の値（奥と書かれた方のバイアス電圧値）
colorplot = 'flow_t';  %【input】カラープロット種類('flow_t')

% ----------------------- %
% 実験ログからショット番号に対応するマッハプローブの位置を取得
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';  %スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,MP.date);
MP.rlist=T.MachProbeRPosition_cm_(MP.shotlist)*10; %マッハプローブr座標[mm]

MPset.K = 2.5;  % 比例定数
k = 1.38e-23;   % [J/K]
q = 1.60e-19;   % [C]
Te = 10 * q / k; % [K]　% T_e：電子温度→青山さんのp.35を参照 5~8[eV]
mi = 1.672e-27; % [kg]
MPset.cs = sqrt(2 * k * Te / mi);

% MachProbeの計算結果を格納
MPdata2D = MP_calc(machpathname, MP, MPset);
%%
% ch = 6; % チャンネル(z位置)の指定
% plot_MP_r_flow(MP, MPdata2D, ch)
% MP_plot(MP,MPdata2D,colorplot, MP.date)

% 磁気面と共にフローをプロット
% Case-O: 磁気面
load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240611047.mat');
% Case-I: 磁気面
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240610018.mat');
MP_plot_with_psi(MP,MPdata2D,colorplot, MP.date, data2D, grid2D)
% movie_MP(machpathname,MP,MPdata2D,colorplot)


%%
% filename = strcat(pathname.pcbprocessed,'\a039_',num2str(shot(1)),'.mat');
% save(filename,"data2D","grid2D","shot");



%%%%% 以下はB4卒論用

% clear all
% clc

% MP.date = 240204;%【input】マッハプローブ計測日
% MP.shotlist = [6:12 14 15 17:35 37:41];% データとして採用するショット番号
% MP.delay = 4;% トリガタイミングの時間ずれ（オシロスコープ）
% MP.start = 469 + MP.delay;%【input】プロット開始時刻[us]
% MP.shot = 5; % 同じ条件で打ったショットの最大値
% MP.scandt = 0.1; % 計算の刻み幅

% % 測定データが入っているフォルダ
% machpathname = 'C:\\Users\\w-har\\OneDrive - The University of Tokyo\\Lab\\pcb_experiment\MachProbe_data';

% MP.mesh = 9;%【input】マッハプローブ補間メッシュ数（z方向）
% MP.meshr = 9;%【input】マッハプローブ補間メッシュ数（r方向）
% MP.trange = 460:MP.scandt:550;%【input】計算時間範囲
% MP.tate = 6;%【input】プロット枚数(縦)
% MP.yoko = 2;%【input】プロット枚数(横)
% MP.dt = 5;%【input】プロット時間間隔[us] plot_MP用
% MP.dt2 = 3;%【input】プロット時間間隔[us] plot_MP_r_flow用
% colorplot = 'flow_t';%【input】カラープロット種類('flow_t')

% % ----------------------- %
% % 実験ログからショット番号に対応するマッハプローブの位置を取得
% DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
% T=getTS6log(DOCID);
% node='date';
% T=searchlog(T,node,MP.date);
% MP.rlist=T.MachProbeRPosition_cm_(MP.shotlist)*10; %マッハプローブr座標[mm]

% MPset.K = 2.5; % 比例定数
% k = 1.38e-23;    % [J/K]
% q = 1.60e-19;  % [C]
% Te = 5 * q / k; % [K]　% T_e：電子温度→青山さんのp.35を参照 5~8[eV]
% mi = 1.672e-27;   % [kg]
% MPset.cs = sqrt(2 * k * Te / mi);


% MPdata2D = cal_MP_ver2(machpathname, MP, MPset);
% %%
% ch = 5; % チャンネル(z位置)の指定
% % plot_MP_r_flow(MP, MPdata2D, ch)
% plot_MP(MP,MPdata2D,colorplot)
% % movie_MP(machpathname,MP,MPdata2D,colorplot)