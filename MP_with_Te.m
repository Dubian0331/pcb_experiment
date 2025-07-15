%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% マッハプローブ
% 流速をプロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
% addpath '/Users/rsomeya/Documents/lab/matlab/common';
% run define_path.m
clc

% データとして採用するショット番号: case-I, 240610
% MP.date = 240610;%【input】マッハプローブ計測日
% MP.shotlist = [18, 20, 23, 26:27, 31:33, 35:38, 40:41, 44, 47];

% データとして採用するショット番号: case-O, 240611
MP.date = 240611;%【input】マッハプローブ計測日
MP.shotlist = [47, 48, 51:52, 54:60, 62, 65:70];

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
MPdata2D = MP_Iratio_calc(machpathname, MP, MPset);

%%

% 1. 補間グリッド(zq)の中で、z=0に最も近い列のインデックスを探す
%    zqは(r,z)のグリッドなので、どの行を使ってもz座標は同じです。ここでは1行目を使います。
[~, z_index] = min(abs(MPdata2D.zq(1, :)));

fprintf('補間グリッド上でz=0に最も近いのは、%d番目の列のデータです。\n', z_index);
fprintf('その位置のZ座標: %f m\n', MPdata2D.zq(1, z_index));

% 2. I_ratio_forplotから、そのインデックスに対応するデータを全て抜き出す
%    I_ratio_forplotの次元は (時間, r座標, z座標) です。
I_ratio_interpolated_at_z0 = MPdata2D.I_ratio_forplot(:, :, z_index);

% これで (時間 x r座標) の2次元行列が得られます。
disp('z=0における補間後のI_ratioデータを抽出しました。');
disp('変数 I_ratio_interpolated_at_z0 に格納されています。');

%%
% --- 1. 電子温度(Te)のMATファイルをロード ---
fprintf('電子温度(Te)のデータをロードしています...\n');
if MP.date == 240610 % Case-I
    te_mat_path = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240610\triple_data2D_Case-I.mat';
elseif MP.date == 240611 % Case-O
    te_mat_path = 'C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\triple_probe\mat\240611\triple_data2D_Case-O.mat';
else
    error('対応するTeの.matファイルが見つかりません。');
end

if exist(te_mat_path, 'file')
    load(te_mat_path, 'triple_data2D', 'R_values', 'time_values');
    % Teデータは(R座標, 時間)の次元で格納されている
    Te_data_3D = triple_data2D.Te_avg;
    Te_data_eV = squeeze(Te_data_3D); % (R, 1, Time) -> (R, Time)
    Te_R_axis = R_values;       % Teの元のR軸
    Te_time_axis = time_values; % Teの元の時間軸
    fprintf('Teデータのロード完了。\n');
else
    error('指定されたTeの.matファイルが見つかりません: %s', te_mat_path);
end


% --- 2. TeのグリッドをI_ratioのグリッドに合わせる (補間) ---
% I_ratioデータのグリッド軸を取得
flow_time_axis = MP.trange;
flow_R_axis = MPdata2D.rq(:, 1); % or r = linspace(min(MP.rlist),max(MP.rlist),MP.meshr)*1E-3;

% 元のTeデータ用のグリッドを作成
[Te_TimeGrid_orig, Te_RGrid_orig] = meshgrid(Te_time_axis, Te_R_axis);

% 補間先のターゲットグリッドを作成
[Flow_TimeGrid_target, Flow_RGrid_target] = meshgrid(flow_time_axis, flow_R_axis);

% interp2を使ってTeデータをI_ratioのグリッドに補間
fprintf('Teデータをフローのグリッドに補間しています...\n');
Te_aligned_eV = interp2(Te_TimeGrid_orig, Te_RGrid_orig, Te_data_eV, Flow_TimeGrid_target, Flow_RGrid_target, 'linear');

% 補間後のTeは(R座標, 時間)なので、(時間, R座標)に転置する
Te_aligned_eV = Te_aligned_eV.';


% --- 3. フロー速度の計算 ---
% 物理定数
k = 1.38e-23;   % [J/K]
mi = 1.672e-27; % [kg] (水素イオン)
K2ev = 11604.5250061657;
K_factor = 2.5; % 比例定数

% 各点でのTe(eV)からイオン音速Cs(m/s)を計算
Te_aligned_eV(Te_aligned_eV < 0) = 0.1;
Te_aligned_K = Te_aligned_eV * K2ev; % eVをKelvinに変換
Cs_2D = sqrt(2 * k * Te_aligned_K / mi);

% 各点でのフロー速度を計算
% I_ratio_interpolated_at_z0 と Cs_2D は同じ (時間, r) の次元を持つ
V_flow_2D = (Cs_2D / K_factor) .* I_ratio_interpolated_at_z0;


% --- 4. フロー速度の2Dカラープロット ---
figure('Name', ['Flow Velocity at z=0 for Shot ' num2str(MP.shotlist(1)) '-' num2str(MP.shotlist(end))], 'NumberTitle', 'on');
[T_grid_plot, R_grid_plot] = meshgrid(flow_time_axis, flow_R_axis);

% contourfは(X, Y, Z)の順で引数をとる。Zは(Yの次元, Xの次元)である必要があるため、V_flowを転置する
contourf(T_grid_plot, R_grid_plot, (V_flow_2D/1000)', 100, 'LineColor', 'none');

% カラーバーとラベル
h = colorbar;
ylabel(h, 'Flow Velocity [km/s]');
clim([-40, 40]); % カラーマップの範囲
colormap('jet');

% 軸ラベルとタイトル
xlabel('Time (us)');
ylabel('R [m]');
% title('Flow Velocity at z=0');
ylim([0.1, 0.3]);
% grid on;

disp('フロー速度のプロットが完了しました。');


%%
% ch = 6; % チャンネル(z位置)の指定
% plot_MP_r_flow(MP, MPdata2D, ch)
% MP_plot(MP,MPdata2D,colorplot, MP.date)

% 磁気面と共にフローをプロット
% Case-O: 磁気面
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240611047.mat');
% Case-I: 磁気面
% load('C:\Users\w-har\OneDrive - The University of Tokyo\Lab\pcb_experiment\processed_data\240610018.mat');
% MP_plot_with_psi(MP,MPdata2D,colorplot, MP.date, data2D, grid2D)
% movie_MP(machpathname,MP,MPdata2D,colorplot)

