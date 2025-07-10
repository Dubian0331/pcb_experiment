%%%%%%%%%%%%%%%%%%%%%%%%
%200ch/280ch用新規pcbプローブと装置の磁場信号チェック
%%%%%%%%%%%%%%%%%%%%%%%%

function [] = probecheck_script(ok_bz,ok_bt,bz,bt)
    figure_switch = ["on","on"];%bz, bt
    row = 14;%プローブ本数＝グラフ出力時の縦に並べる個数(280chの場合は14、200chの場合は10）
    col = 10;%グラフ出力時の横に並べる個数(r方向10ch per probe)
    upper_lim_bz = 0.05;%縦軸プロット領域（b_z上限）
    lower_lim_bz = -0.05;%縦軸プロット領域（b_z下限）
    upper_lim_bt = 0.2;%縦軸プロット領域（b_t上限）
    lower_lim_bt = -0.2;%縦軸プロット領域（b_t下限）
    t_start=400;%横軸プロット領域（開始時間）
    t_end=500;%横軸プロット領域（終了時間）
    
    f1=figure(Visible=figure_switch(1));
    f1.WindowState = 'maximized';
    for i=1:row
        for j=1:col
            subplot(row,col,(i-1)*col+j)
            if ok_bz(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
                plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j));
            else %NGなチャンネルは赤色点線でプロット
                plot(t_start:t_end,bz(t_start:t_end,col*(i-1)+j),'r:')
            end   
            title(num2str(2.*(col*(i-1)+j)-1));
            xticks([t_start t_end]);
            ylim([lower_lim_bz upper_lim_bz]);
        end
    end
    sgtitle(['Bz signal probe1-',num2str(row)]) 
    
    f2=figure(Visible=figure_switch(2));
    f2.WindowState = 'maximized';
    for i=1:row
        for j=1:col
            subplot(row,col,(i-1)*col+j)
            if ok_bt(col*(i-1)+j)==1 %okなチャンネルはそのままプロット
                plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j))
            else %NGなチャンネルは赤色点線でプロット
                plot(t_start:t_end,bt(t_start:t_end,col*(i-1)+j),'r:')
            end   
            title(num2str(2.*(col*(i-1)+j)));
            xticks([t_start t_end]);
            ylim([lower_lim_bt upper_lim_bt]);
        end
    end
    sgtitle(['Bt signal probe1-',num2str(row)]) 
    
    hidden = find(figure_switch == 'off');
    figures = [f1,f2];
    for i = hidden
        close(figures(i));
    end

end