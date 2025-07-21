% 実験ログからショット番号に対応するマッハプローブの位置を取得
date = 240610;
shotlist = [1:48];
DOCID='1wG5fBaiQ7-jOzOI-2pkPAeV6SDiHc_LrOdcbWlvhHBw';%スプレッドシートのID
T=getTS6log(DOCID);
node='date';
T=searchlog(T,node,date);
rlist=T.tripleProbeRPosition_cm_____________________0_5cm(shotlist)*0.01; %マッハプローブr座標[m]
a039_lists = T.a039(shotlist);

simple_list = cat(2, a039_lists, rlist);