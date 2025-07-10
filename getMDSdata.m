function [x]=getMDSdata(dtacq_num,shot,tfshot)
%%%指定したデジタイザの192ch分のデータを読み込む＋オフセット＋TF差し引きをする関数
%%%【input】dtacq:38/39/40,
%%%shot:dtacqのshot番号,tfshot:TFoffsetに対応するdtacqのshot番号、ない場合は0
clear x y
if dtacq_num==38
    ch_num=128;
else
    ch_num=192;
end
post=1000;%t=0からの計測時間[us]
dtacq=strcat('a',num2str(dtacq_num,'%03i'));%a038などの形式の文字列へ変換
x=zeros(post,ch_num);
y=x;

import MDSplus.*
%ツリーのdatafileがあるフォルダのパスをtreename_pathという形で環境変数に設定(a038_path, a039_path, a040_path)
%mdsipのポートに接続して各デジタイザののツリーを開く。
mdsconnect('192.168.1.140');
mdsopen(dtacq, shot); 

for i=1:ch_num
    %各チャンネルにおいて「.AI:CHXXX」というノードを指定するためのノード名を作る
    chname=".AI:CH"+num2str(transpose(i),'%03i');
    % num2strで数値データをstrデータに変換。この時'%03i'で左側を(0)で埋めた(3)桁の整数(i)という形を指定できる。
    x(:,i)=mdsvalue(chname);
    %データがとれていないときエラーメッセージが多分237文字で帰ってくるので、1000以下の要素はデータなしとしてリターンする
    if numel(x(:,i)) <1000
        return
    end
    x(:,i)=x(:,i)-x(1,i);% オフセット調整
end
if tfshot>0 
    mdsopen(dtacq, tfshot);
    for i=1:ch_num
    %各チャンネルにおいて「.AI:CHXXX」というノードを指定するためのノード名を作る
    chname=".AI:CH"+num2str(transpose(i),'%03i');
    y(:,i)=mdsvalue(chname);
    if numel(y(:,i)) <1000
        y=zeros(post,ch_num);
    end
    y(:,i)=y(:,i)-y(1,i);% オフセット調整
    end
    x=x-y;% TFノイズを差し引いたもの
end

end