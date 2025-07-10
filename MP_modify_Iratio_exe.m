% --- I_ratio修正版を一括生成する処理 ---

date = 240610;
src_dir = 'C:/Users/w-har/OneDrive - The University of Tokyo/Lab/pcb_experiment/MachProbe_data';
dst_dir = 'C:/Users/w-har/OneDrive - The University of Tokyo/Lab/pcb_experiment/MachProbe_data';

for shotnum = 35:47
    try
        MP_modify_Iratio(date, shotnum, src_dir, dst_dir);
    catch ME
        disp(['修正版作成でエラー: shotnum = ', num2str(shotnum)]);
        disp(ME.message);
        continue;
    end
end