function [x]=getvalue(shot,tfshot)
%%%128ch分のデータを読み込む＋オフセット＋TF差し引きをする関数
ch_num=128;

import MDSplus.*
%ツリーのdatafileがあるフォルダのパスをtreename_pathという形で環境変数に設定
%setenv('a038_path','F:\a038');
setenv('a038_path','K:\mnt\fourier\a038');
%mdsipのポートに接続してa038のツリーを開く。
mdsconnect('192.168.1.140');
mdsopen('a038', shot); %a038ツリーの1768ショットを開く
%「.AI:CH001」というノードを指定してxにCH001のデータを出力してプロット


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
if exist("tfshot")==1 
    mdsopen('a038', tfshot);
    for i=1:ch_num
    %各チャンネルにおいて「.AI:CHXXX」というノードを指定するためのノード名を作る
    chname=".AI:CH"+num2str(transpose(i),'%03i');
    y(:,i)=mdsvalue(chname);
    y(:,i)=y(:,i)-y(1,i);% オフセット調整
    end
    x=x-y;% TFノイズを差し引いたもの
end


end

