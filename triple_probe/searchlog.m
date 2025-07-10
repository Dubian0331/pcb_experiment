%%logから日付や語句検索でショットを探す関数
% node='テーブルのヘッダー'：　検索するもの　'operator','date',など
% pat=数値 or str :　抽出するパターン
function output=searchlog(inputtable,node,pat)
%%名前,commentで検索(~を含むの検索)
if ischar(pat)
    rows = contains(inputtable.(node),pat,'IgnoreCase',true);
    %%日付,shotなどで検索(一致検索)
else
    rows = inputtable.(node) == pat ;
end
%%dateが空白のshotも取り込む
if node == 'date'
    max_row = find(rows,1,'last');%dateが日付と一致する最終行
    if max_row < size(rows,1)
        while isnan(inputtable.date(max_row+1)) && isnumeric(inputtable.shot(max_row+1))
            rows(max_row+1) = 1;
            max_row = max_row + 1;
            if max_row >= size(rows,1)
                break
            end
        end
    end
end
output = inputtable(rows,:);
end