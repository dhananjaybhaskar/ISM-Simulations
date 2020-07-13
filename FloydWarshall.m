function [D, next] = FloydWarshall(D)

    next = zeros(size(D));
    len = length(D);
    
    for i = 1:len
        for j = 1:len
            if(D(i,j) ~= inf)
                next(i, j) = j;
            end
        end
    end
    
    for k = 1:len
        for i = 1:len
            for j = 1:len
                dist = D(i, k) + D(k, j);
                if(D(i, j) > dist && (D(i,k) ~= 0 && D(k,j) ~= 0))
                    D(i, j) = dist;
                    next(i,j) = next(i, k);
                end
            end
        end
    end
    
end
