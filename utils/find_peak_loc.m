function [peak_loc] = find_peak_loc(signal, loc, L)
idx_start = 1;
idx_end = 1;


loc(loc<L) = [];
loc(loc>length(signal)-L)= [];

peak_loc1 = [];
while(idx_end<length(loc))
    while loc(idx_end+1)-loc(idx_start)<L 
        idx_end = idx_end+1;
        if idx_end == length(loc)
            break;
        end
    end
    
    temp = idx_start:idx_end;

    temp_signal = signal(loc(temp));
    [~, temp_idx] = max(temp_signal);
    if size(temp_idx)<L
    peak_loc1 = [peak_loc1, loc(temp_idx+idx_start-1)];
    end
  
    idx_start = idx_end+1;
    idx_end = idx_start;
end

peak_loc = nan(1,length(peak_loc1));
for n = 1:length(peak_loc1)
    temp = signal(peak_loc1(n)-ceil(L/2)+1:peak_loc1(n)+ceil(L/2));
    [~, temp_idx] = max(temp);
    peak_loc(n) = temp_idx+peak_loc1(n)-ceil(L/2);
end
peak_loc = unique(peak_loc);
end
