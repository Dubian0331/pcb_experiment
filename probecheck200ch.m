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

pathname.rawdata='/Users/uebokoki/Library/CloudStorage/GoogleDrive-7162323547@edu.k.u-tokyo.ac.jp/My Drive/lab/lab_data/raw_data';%dtacqのrawdataの保管場所

%%%%実験オペレーションの取得
%直接入力の場合
dtacqlist=39;
shotlist=3640;%【input】dtacqの保存番号
tfshotlist=3638;
date = 240116;%【input】計測日
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
filename=strcat(pathname.rawdata,'/rawdata_dtacq',num2str(dtacq_num),'_shot',num2str(shot),'_tfshot',num2str(tfshot),'.mat');
% filename=strcat(pathname.rawdata,'rawdata_noTF_dtacq',num2str(d_tacq),'.mat');
load(filename,'rawdata');%1000×192

%正しくデータ取得できていない場合はreturn
if numel(rawdata)< 500
    return
end

%較正係数のバージョンを日付で判別
sheets = sheetnames('coeff200ch.xlsx');
sheets = str2double(sheets);
sheet_date=max(sheets(sheets<=date));

C = readmatrix('coeff200ch.xlsx','Sheet',num2str(sheet_date));
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

b=rawdata.*coeff';%較正係数RC/NS
b=b.*P';%極性揃え
b=smoothdata(b,1);

%デジタイザchからプローブ通し番号順への変換
bz=zeros(1000,100);
bt=bz;
ok_bz=true(1,100);
ok_bt=ok_bz;

windowSize = 3;
bb = (1/windowSize)*ones(1,windowSize);
aa = 1;
for i=1:192
    b(:,i) = filter(bb,aa,b(:,i));
    if rem(ch(i),2)==1
        bz(:,ceil(ch(i)/2))=b(:,i);
        ok_bz(ceil(ch(i)/2))=ok(i);

    elseif rem(ch(i),2)==0
        bt(:,ceil(ch(i)/2))=b(:,i);
        ok_bt(ceil(ch(i)/2))=ok(i);
    end
end
ok_bz_plot=ok_bz;
ok_bt([4 5 6 7 8 9 10 15 21 27 30 42 43 49 53 69 84 87 92 94 95 96 97 98 99 100]) = false;

% 221219ver
% Pcheck=[1	-1	1	1	-1	-1	-1	1	1	1	1	1	1	-1	-1	0	-1	-1	1	-1	1	0	1	-1	-1	1	-1	-1	1	-1	-1	1	-1	1	-1	1	-1	-1	-1	1	-1	-1	1	1	-1	-1	-1	-1	-1	-1	1	-1	-1	-1	-1	-1	-1	-1	-1	-1	1	1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	1	-1	1	-1	-1	-1	-1	1	-1	-1	1	-1	-1	1	-1	-1	1	1	1	-1	1	-1	-1	1	1	1	-1];
% bz=bz.*Pcheck;
% bz(:,[37 47 57])=-bz(:,[37 47 57]);
% bz(:,70)=-bz(:,70);
% for i=0:9
%     bz(:,[7 8 9 10]+10.*i)=-bz(:,[7 8 9 10]+10.*i);
% end
% ok_bz([5 6 11 16 20 22 31 33 39 49 63 66 71 72 79 80 95 100])=false;
% ok_bz(49)=true;
% bz(:,49)=-bz(:,49);

%生信号描画用パラメータ
r = 5;%プローブ本数＝グラフ出力時の縦に並べる個数
col = 10;%グラフ出力時の横に並べる個数
y_upper_lim = 0.1;%3e-3;%0.1;%縦軸プロット領域（b_z上限）
y_lower_lim = -0.1;%3e-3;%-0.1;%縦軸プロット領域（b_z下限）
t_start=400;%430;%455;%横軸プロット領域（開始時間）
t_end=600;%550;%横軸プロット領域（終了時間）
% r_ch=col1+col2;%r方向から挿入した各プローブのチャンネル数

f1=figure;
f1.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz_plot(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe1-5')

f2=figure;
f2.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bz_plot(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bz(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)-1));
        xticks([t_start t_end]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bz signal probe6-10')

f3=figure;
f3.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i-1)+j)));
        xticks([t_start t_end]);
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bt signal probe1-5')

f4=figure;
f4.WindowState = 'maximized';
for i=1:r
    for j=1:col
        subplot(r,col,(i-1)*col+j)
        if ok_bt(col*(i+r-1)+j)==1 %okなチャンネルはそのままプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j))
        else %NGなチャンネルは赤色点線でプロット
            plot(t_start:t_end,bt(t_start:t_end,col*(i+r-1)+j),'r:')
        end   
        title(num2str(2.*(col*(i+r-1)+j)));
        xticks([t_start t_end]);
        %ylim([-0.2 0.2]);
        ylim([y_lower_lim y_upper_lim]);
    end
end
sgtitle('Bt signal probe6-10')

% saveas(gcf,strcat(pathname.save,'\date',num2str(date),'_dtacq',num2str(d_tacq),'_02','.png'))
% close

%横軸z, 縦軸Bzのプロット
% f5=figure;
% f5.WindowState = 'maximized';
% t=465;
% for i=1:10
%     zline=(1:10:91)+(i-1);
%     bz_zline=bz(t,zline);
%     bz_zline(ok_bz(zline)==false)=NaN;
%     plot([-0.17 -0.1275 -0.0850 -0.0315 -0.0105 0.0105 0.0315 0.0850 0.1275 0.17],bz_zline,'-*')
%     clear bz_zline
%     hold on
% end
% hold off
% xlabel('z [m]')
% ylabel('Bz')
% yline(0,'k--')
% title(strcat('t=',num2str(t),' us'))
% legend('r1','r2','r3','r4','r5','r6','r7','r8','r9','r10',Location='eastoutside')

% bz1=[bz(t,1) bz(t,11) bz(t,21) bz(t,31) bz(t,41) bz(t,51) bz(t,61) bz(t,71) bz(t,81) bz(t,91)];
% bz2=[bz(t,2) bz(t,12) bz(t,22) bz(t,32) bz(t,42) bz(t,52) bz(t,62) bz(t,72) bz(t,82) bz(t,92)];
% bz3=[bz(t,3) bz(t,13) bz(t,23) bz(t,33) bz(t,43) bz(t,53) bz(t,63) bz(t,73) bz(t,83) bz(t,93)];
% bz4=[bz(t,4) bz(t,14) bz(t,24) bz(t,34) bz(t,44) bz(t,54) bz(t,64) bz(t,74) bz(t,84) bz(t,94)];
% bz5=[bz(t,5) bz(t,15) bz(t,25) bz(t,35) bz(t,45) bz(t,55) bz(t,65) bz(t,75) bz(t,85) bz(t,95)];
% bz6=[bz(t,6) bz(t,16) bz(t,26) bz(t,36) bz(t,46) bz(t,56) bz(t,66) bz(t,76) bz(t,86) bz(t,96)];
% bz7=[bz(t,7) bz(t,17) bz(t,27) bz(t,37) bz(t,47) bz(t,57) bz(t,67) bz(t,77) bz(t,87) bz(t,97)];
% bz8=[bz(t,8) bz(t,18) bz(t,28) bz(t,38) bz(t,48) bz(t,58) bz(t,68) bz(t,78) bz(t,88) bz(t,98)];
% bz9=[bz(t,9) bz(t,19) bz(t,29) bz(t,39) bz(t,49) bz(t,59) bz(t,69) bz(t,79) bz(t,89) bz(t,99)];
% bz10=[bz(t,10) bz(t,20) bz(t,30) bz(t,40) bz(t,50) bz(t,60) bz(t,70) bz(t,80) bz(t,90) bz(t,100)];

end
