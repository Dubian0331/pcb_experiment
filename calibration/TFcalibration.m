%%%プローブをBt方向計測に90度回転させ、a038とa039の感度の相対値を較正
%TF onlyの場合のBtは1/rに比例→同じr位置での値の比較
%2023/01/10計測：a039のch5, 9は回転させにくかったのでBz計測方向のまま

pathname.ts3u='C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\rgwdata';%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先

pathname.rawdata38=getenv('rawdata038_path');%dtacq a038のrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所

pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所

dtacq_num=[38 39];
shot=[10367 308];%【input】dtacqの保存番号
tfshot=[0 0];
date = 230110;%【input】計測日
rgwshot=2;

filename1=strcat(pathname.rawdata,'\rawdata_dtacq',num2str(dtacq_num(1)),'_shot',num2str(shot(1)),'_tfshot',num2str(tfshot(1)),'.mat');
load(filename1,'rawdata');%a038
rawdata1=rawdata;
clear rawdata
filename2=strcat(pathname.rawdata,'\rawdata_dtacq',num2str(dtacq_num(2)),'_shot',num2str(shot(2)),'_tfshot',num2str(tfshot(2)),'.mat');
load(filename2,'rawdata');%a039
rawdata2=rawdata;
clear rawdata

%正しくデータ取得できていない場合はreturn
if numel(rawdata1)< 500||numel(rawdata2)<500
    return
end

%較正係数のバージョンを日付で判別
sheets1 = sheetnames('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff125ch.xlsx');
sheets1 = str2double(sheets1);
sheet_date1=max(sheets1(sheets1<=date));
sheets2 = sheetnames('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff200ch.xlsx');
sheets2 = str2double(sheets2);
sheet_date2=max(sheets2(sheets2<=date));
C1 = readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff125ch.xlsx','Sheet',num2str(sheet_date1));
C2 = readmatrix('C:\Users\kuru1\OneDrive - g.ecc.u-tokyo.ac.jp\labo\experiment\coeff200ch.xlsx','Sheet',num2str(sheet_date2));

%a038
ok1 = logical(C1(:,14));
P1=C1(:,13);
coeff1=C1(:,12);
zpos1=C1(:,9);
rpos1=C1(:,10);
probe_num1=C1(:,5);
probe_ch1=C1(:,6);
ch1=C1(:,7);

b1=rawdata1.*coeff1';%較正係数RC/NS
b1=b1.*P1';%極性揃え
b1=double(b1);

%デジタイザchからプローブ通し番号順への変換
bz1=zeros(1000,126);
ok_bz1=true(100,1);
rpos_bz1=zeros(126,1);
zpos_bz1=rpos_bz1;

for i=1:128
    if ch1(i)>0
        bz1(:,ch1(i))=b1(:,i);
        ok_bz1(ch1(i))=ok1(i);
        rpos_bz1(ch1(i))=rpos1(i);
        zpos_bz1(ch1(i))=zpos1(i);
    end
end
bz1(:,63)=[];
ok_bz1(63)=[];
rpos_bz1(63)=[];
zpos_bz1(63)=[];

%a039
ok2 = logical(C2(:,14));
P2=C2(:,13);
coeff2=C2(:,12);
zpos2=C2(:,9);
rpos2=C2(:,10);
probe_num2=C2(:,5);
probe_ch2=C2(:,6);
ch2=C2(:,7);

b2=rawdata2.*coeff2';%較正係数RC/NS
b2=b2.*P2';%極性揃え

%デジタイザchからプローブ通し番号順への変換
bz2=zeros(1000,100);
bt2=bz2;
ok_bz2=false(100,1);
ok_bt2=ok_bz2;
zpos_bz2=zeros(100,1);
rpos_bz2=zpos_bz2;
zpos_bt2=zpos_bz2;
rpos_bt2=zpos_bz2;

for i=1:192
    if rem(ch2(i),2)==1
        bz2(:,ceil(ch2(i)/2))=b2(:,i);
        ok_bz2(ceil(ch2(i)/2))=ok2(i);
        zpos_bz2(ceil(ch2(i)/2))=zpos2(i);
        rpos_bz2(ceil(ch2(i)/2))=rpos2(i);
    elseif rem(ch2(i),2)==0
        bt2(:,ch2(i)/2)=b2(:,i);
        ok_bt2(ceil(ch2(i)/2))=ok2(i);
        zpos_bt2(ceil(ch2(i)/2))=zpos2(i);
        rpos_bt2(ceil(ch2(i)/2))=rpos2(i);
    end
end

%立ち上げノイズ削減
for i=1:100
    bz2(:,i)=bz2(:,i)-mean(bz2(10:30,i));
end
for i=1:125
    bz1(:,i)=bz1(:,i)-mean(bz1(10:70,i));
end

t=1:1000;

ok_bz2([41:50 81:90 35 36 27])=false;%Bt計測に回していないプローブch＋死んでいるch追加
r_2=[1 3 4 5 6 7];%r=0.06,0.12,0.15,0.18,0.21,0.24
r_1=[1 4 6 12 18 24];
bz_max=zeros(6,2);

for i=1:6
    ok2_cal=ok_bz2(r_2(i)+(0:10:90));
    bz2_r=bz2(:,r_2(i)+(0:10:90));
    bz2_r=bz2_r(:,ok2_cal);
    bz_max(i,2)=mean(min(bz2_r(100:end,:),[],1));
    ok1_cal=ok_bz1(r_1(i)+(0:25:100));
    bz1_r=bz1(:,r_1(i)+(0:25:100));
    bz1_r=bz1_r(:,ok1_cal);
    bz_max(i,1)=mean(min(bz1_r(100:end,:),[],1));
end

bz_coeff=bz_max(:,1)./bz_max(:,2);
bz_coeff_mean=mean(bz_coeff);
%a038はa039の1.1976倍


for k=1:10
    if ok_bz2(r_ch+(k-1).*10)==1
        mean(:,bz2(:,r_ch+(k-1).*10))
    end
end

figure
for i=1:2
    for j=1:5
        r_ch=(i-1).*5+j;
        subplot(2,5,r_ch)
        for k=1:10
            if ok_bz2(r_ch+(k-1).*10)==1
                plot(t,bz2(:,r_ch+(k-1).*10))
                hold on
            end
        end
        yline(0,'k--')
        hold off
        ylim([-0.5 0.1])
        title(num2str(r_ch))
    end
end

ok_bz1([63 96 119 121 44 107])=false;%Bt計測に回していないプローブch＋死んでいるch追加
figure
for i=1:2
    for j=1:6
        r_ch=(i-1).*6+j;
        subplot(2,6,r_ch)
        for k=1:5
            if ok_bz1(r_ch+(k-1).*25)==1
                plot(t,bz1(:,r_ch+(k-1).*25))
                hold on
            end
        end
        yline(0,'k--')
        hold off
        ylim([-0.5 0.1])
        title(num2str(r_ch))
    end
end
figure
for i=1:2
    for j=1:7
        r_ch=(i-1).*7+j+12;
        if r_ch==26
            return
        else
        subplot(2,7,r_ch-12)
        for k=1:5
            if ok_bz1(r_ch+(k-1).*25)==1
                plot(t,bz1(:,r_ch+(k-1).*25))
                hold on
            end
        end
        yline(0,'k--')
        hold off
        ylim([-0.5 0.1])
        title(num2str(r_ch))
        end
    end
end


% figure
% plot(t,bz2(:,1),t,bz2(:,11),t,bz2(:,21),t,bz2(:,51),t,bz2(:,91))
% yline(0,'k--')
% ylim([-0.5 0.1])

% % rgwデータの読み込み->coil current plot用
% filepath.rgw=strcat(pathname.ts3u, '/', string(date),'/' ...
%     ,string(date),num2str(rgwshot,'%03i'),'.rgw');
% rgwdata = importfile(filepath.rgw);
% % rgwdata.ch6=rgwdata.ch6*(65.7146)*12;
% % rgwdata.ch7=rgwdata.ch7*(-57.6030)*3;
% % rgwdata.ch8=rgwdata.ch8*(-187.170)*3-50;
% %各ショットのTF,PF1,PF2のcoil current
% % figure
% % plot(repmat([0:0.1:0.1*10000],4,1)',rgwdata{2:end,:})
% % legend('ch1','ch2','ch5','ch7')
% % ylabel('coil current [kA]')
% % % ylim([-200 150])
% % xlim([0 1000])
% % title(strcat('IDX',num2str(IDX),': Cr',num2str(T.Crowbar_us_(IDX))))
% 
% function rgwdata = importfile(filename, dataLines)
% %IMPORTFILE テキスト ファイルからデータをインポート
% %  UNTITLED = IMPORTFILE(FILENAME) は既定の選択に関してテキスト ファイル FILENAME
% %  からデータを読み取ります。  データを table として返します。
% %
% %  UNTITLED = IMPORTFILE(FILE, DATALINES) はテキスト ファイル FILENAME
% %  の指定された行区間のデータを読み取ります。DATALINES
% %  を正の整数スカラーとして指定するか、行区間が不連続の場合は正の整数スカラーからなる N 行 2 列の配列として指定します。
% %
% %  例:
% %  Untitled = importfile("X:\results\ts-3u\211222\211222006.rgw", [1, Inf]);
% %
% %  READTABLE も参照してください。
% %
% % MATLAB からの自動生成日: 2022/01/15 13:57:51
% 
% %% 入力の取り扱い
% % dataLines が指定されていない場合、既定値を定義します
% if nargin < 2
%     dataLines = [1, Inf];
% end
% 
% %% インポート オプションの設定およびデータのインポート
% opts = delimitedTextImportOptions("NumVariables", 18);
% 
% % 範囲と区切り記号の指定
% opts.DataLines = dataLines;
% opts.Delimiter = "\t";
% 
% % % 列名と型の指定
% % opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "ch6", "ch7", "ch8", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"];
% % opts.SelectedVariableNames = ["ch6", "ch7", "ch8"];
% % opts.VariableTypes = ["string", "string", "string", "string", "string", "string", "string", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string"];
% % 列名と型の指定
% opts.VariableNames = ["Var0", "ch1", "ch2", "Var3", "Var4", "ch5", "Var6", "ch7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17"];
% opts.SelectedVariableNames = ["ch1", "ch2", "ch5", "ch7"];
% opts.VariableTypes = ["string", "double", "double", "string", "string", "double", "string", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];
% 
% % ファイル レベルのプロパティを指定
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% 
% % % 変数プロパティを指定
% % opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"], "WhitespaceRule", "preserve");
% % opts = setvaropts(opts, ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"], "EmptyFieldRule", "auto");
% % 変数プロパティを指定
% opts = setvaropts(opts, ["Var0", "Var3", "Var4", "Var6", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17"], "WhitespaceRule", "preserve");
% opts = setvaropts(opts, ["Var0", "Var3", "Var4", "Var6", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17"], "EmptyFieldRule", "auto");
% 
% 
% % データのインポート
% rgwdata = readtable(filename, opts);
% 
% end



