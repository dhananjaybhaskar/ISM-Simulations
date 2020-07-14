function [translated_path_idxes] = compute_path(next, i, j, num_cols)

    Path = zeros(2*length(next), 1);
    Path(1) = i;
    idx = 2;
    cur_i = i;
    
    while(cur_i ~= j)
        cur_i = next(cur_i, j);
        Path(i,j,idx) = cur_i;
        idx = idx + 1;
    end
    
    Path = Path(Path ~= 0); 
    Path = Path(:);
    translated_path_idxes = zeros(numel(Path), 2);
    
    for i = 1:numel(Path)
        assert(Path(i) ~= 0);
        translated_path_idxes(i, 2) = mod(Path(i) - 1, num_cols) + 1;
        translated_path_idxes(i, 1) = (Path(i) - translated_path_idxes(i,2))/num_cols + 1;
    end
    
end