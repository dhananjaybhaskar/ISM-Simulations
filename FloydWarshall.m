% function [D,Path, K] = FloydWarshall(D)
%     prevD = D;
%     P = zeros(size(D));
%     for k = 1:length(D)
%       D = min(D,D(:,k) + D(k,:));
%       P(D<prevD) = k;
%       prevD = D;
%     end
%     Path = zeros([size(D) 2*length(D)]);
%     for i = 1:length(D)
%         for j = 1:length(D)
%             cur_j = j;
%             idx = 2*length(D);
%             while(cur_j ~= i && cur_j ~= 0)
%                 Path(i,j,idx) = cur_j;
%                 idx = idx - 1;
%                 cur_j = P(i,cur_j);
%             end
%             Path(i,j,idx) = i;
%         end
%     end
%     K = P;
% end

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
%     Path = zeros(len, len, 2*len);
%     for i = 1:len
%         Path(i,:,1) = i;
%         for j = 1:len
%             idx = 2;
%             cur_i = i;
%             while(cur_i ~= j)
%                 cur_i = next(cur_i, j);
%                 Path(i,j,idx) = cur_i;
%                 idx = idx + 1;
%             end
%         end
%     end
end
