%%%%%%%%%%%%%%%%%%%%%%%%
%200ch用新規pcbプローブと装置の磁場信号極性チェック
%dtacqのshot番号を直接指定する場合
%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先

pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所

pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所

%%%%実験オペレーションの取得
%直接入力の場合
dtacqlist=38;
shotlist=10531;%【input】dtacqの保存番号
tfshotlist=10530;
date = 230119;%【input】計測日
n=numel(shotlist);%計測データ数

% %磁気面出す場合は適切な値を入力、磁場信号のみプロットする場合は変更不要
% d_tacqTF = '';%【input】TFoffsetのdtacq保存番号
% i_EF = 0;%【input】EF電流
% trange=460:490;%【input】計算時間範囲
% n=50; %【input】rz方向のメッシュ数

for i=1:n
    dtacq_num=dtacqlist(i);
    shot=shotlist(i);
    tfshot=tfshotlist(i);
    check_signal(date, dtacq_num, shot, tfshot, pathname); 
end

%%%%%%%%%%%%%%%%%%%%%%%%
%以下、local関数
%%%%%%%%%%%%%%%%%%%%%%%%

function check_signal(date, dtacq_num, shot, tfshot, pathname)
filename=strcat(pathname.rawdata,'\rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
% filename=strcat(pathname.rawdata,'rawdata_noTF_dtacq',num2str(d_tacq),'.mat');
load(filename,'rawdata');%1000×192

%正しくデータ取得できていない場合はreturn
if numel(rawdata)< 500
    return
end

%較正係数のバージョンを日付で判別
sheets = sheetnames('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff125ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));

C = readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff125ch.xlsx','Sheet',num2str(sheet_date));
ok = logical(C(:,14));
P=C(:,13);
coeff=C(:,12);
zpos=C(:,9);
rpos=C(:,10);
probe_num=C(:,5);
probe_ch=C(:,6);
ch=C(:,7);
d2p=C(:,15);
d2bz=C(:,16);
d2bt=C(:,17);

p_ch= readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff125ch.xlsx','Sheet','p_ch');

b=rawdata.*coeff';%較正係数RC/NS
b=b.*P';%極性揃え
b=smoothdata(b,1);

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,126);
ok_bz=true(1,100);

for i=1:128
    if ch(i)>0
        bz(:,ch(i))=b(:,i);
        ok_bz(ch(i))=ok(i);
    end
end
bz(:,63)=[];
ok_bz(63)=[];

% [bz, ok_bz, ok_bz_plot] = ng_replace(bz, ok_bz, sheet_date);


bz_s=bz;
for i=1:125
    bz_s(:,i)=lowpass(bz(:,i),0.4e5,1e6);
end
% bz_s=filloutliers(bz_s,"previous","movmean",5);
ok_bz([5 54 31 8 61 37 62 13 38 16 40 20 45 47 23 48 25])=false;

%生信号描画用パラメータ
r = 5;%プローブ本数＝グラフ出力時の縦に並べる個数
col1 = 12;%1枚目のグラフ出力時の横に並べる個数
col2 = 13;%2枚目のグラフ出力時の横に並べる個数
y_upper_lim = 0.05;%3e-3;%0.1;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.05;%3e-3;%-0.1;%縦軸プロット領域（b_z下限）
t_start=350;%430;%455;%横軸プロット領域（開始時間）
t_end=530;%550;%横軸プロット領域（終了時間）
r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数

f=figure;
f.WindowState = 'maximized';
for i=1:r
    for j=1:col1
        subplot(r,col1,(i-1)*col1+j)
        if ok_bz(r_ch*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
            hold on
            plot(t_start:t_end,bz_s(t_start:t_end,r_ch*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
        end   
        title(num2str(p_ch(i,j)));
        %xlim([t_start t_end]);
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
        %ylim([-0.02 0.04]);
    end
end
% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_01','.png'))
% close

f2=figure;
f2.WindowState = 'maximized';
for i=1:r
    for j=col1+1:col1+col2
        subplot(r,col2,(i-1)*col2+j-col1)
        if ok_bz(r_ch*(i-1)+j)==1 
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j))
            hold on
            plot(t_start:t_end,bz_s(t_start:t_end,r_ch*(i-1)+j))
        else 
            plot(t_start:t_end,bz(t_start:t_end,r_ch*(i-1)+j),'r:')
        end   
        title(num2str(p_ch(i,j)));
        %xlim([t_start t_end]);
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
        %ylim([-0.02 0.04]);
    end
end
% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_02','.png'))
% close

%横軸z, 縦軸Bzのプロット
f5=figure;
f5.WindowState = 'maximized';
t=470;
subplot(3,1,1)
for i=1:8
    zline=(1:25:101)+(i-1);
    bz_zline=bz(t,zline);
    bz_zline(ok_bz(zline)==false)=NaN;
    plot(1:5,bz_zline,'-*')
    clear bz_zline
    hold on
end
hold off
xlabel('z [m]')
ylabel('Bz')
yline(0,'k--')
legend('r1','r2','r3','r4','r5','r6','r7','r8',Location='eastoutside')

subplot(3,1,2)
for i=9:16
    zline=(1:25:101)+(i-1);
    bz_zline=bz(t,zline);
    bz_zline(ok_bz(zline)==false)=NaN;
    plot(1:5,bz_zline,'-*')
    clear bz_zline
    hold on
end
hold off
xlabel('z [m]')
ylabel('Bz')
yline(0,'k--')
legend('r9','r10','r11','r12','r13','r14','r15','r16',Location='eastoutside')

subplot(3,1,3)
for i=17:25
    zline=(1:25:101)+(i-1);
    bz_zline=bz(t,zline);
    bz_zline(ok_bz(zline)==false)=NaN;
    plot(1:5,bz_zline,'-*')
    clear bz_zline
    hold on
end
hold off
xlabel('z [m]')
ylabel('Bz')
yline(0,'k--')
legend('r17','r18','r19','r20','r21','r22','r23','r24','r25',Location='eastoutside')
sgtitle(strcat('t=',num2str(t),' us'))
end
