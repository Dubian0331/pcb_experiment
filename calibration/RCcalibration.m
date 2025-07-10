%%%デジタイザを用いた積分器RC較正データの取得

%%%%%ここが各PCのパス
%【※コードを使用する前に】環境変数を設定しておくか、matlab内のコマンドからsetenv('パス名','アドレス')で指定してから動かす
pathname.ts3u=getenv('ts3u_path');%old-koalaのts-3uまでのパス（mrdなど）
pathname.fourier=getenv('fourier_path');%fourierのmd0（データックのショットが入ってる）までのpath
pathname.NIFS=getenv('NIFS_path');%resultsまでのpath（ドップラー、SXR）
pathname.save=getenv('savedata_path');%outputデータ保存先

pathname.rawdata=getenv('rawdata_path');%dtacqのrawdataの保管場所
pathname.woTFdata=getenv('woTFdata_path');%rawdata（TFoffset引いた）の保管場所
% 
% %計測chなど読み込み
% S=readmatrix("221216RC.xlsx");
% fail=S(:,5);
% S(fail==1,:)=[];
% 
% dtacq_list=S(:,1);
% shot_list=S(:,2);
% int_list=S(:,3);
% fg_list=S(:,4);

%手動入力
dtacq_list=39;
shot_list=9;
int_list=1;
fg_list=128;

n=numel(shot_list);%計測データ数
data=zeros(n,3);%較正係数保管

for i=1:n
    dtacq_num=dtacq_list(i);
    shot=shot_list(i);
    int_ch=int_list(i);
    fg_ch=fg_list(i);
    
    x=getMDSdata(dtacq_num,shot,0);
    int=x(:,int_ch);
    fg=x(:,fg_ch);
    int=int-mean(int);%offsetを0に揃える
    fg=fg-mean(fg);
    RC=rms(int)./rms(fg);

    data(i,1)=dtacq_num;
    data(i,2)=int_ch;
    data(i,3)=RC;

    %信号図示
    t=1:1000;
    figure
    plot(t,int,'r')
    hold on
    plot(t,fg,'b')
    legend('Integrator','Function Generator')
    ylim([-1 1])
    xlabel('t [us]')
    ylabel('V [V]')
    ha1 = gca;
    ha1.LineWidth = 1;
    ha1.FontSize=13;

end

% save(strcat(pathname.save,'\','221216RC.mat'),'data')

