function selectedBands = computeSelectedBands(SimilarityMatrix, TotalBands, NumSubspaces)
    % SimilarityMatrix: Similarity matrix
    % TotalBands: Total number of bands
    % NumSubspaces: Number of subspaces

    selectedBands = linspace(1, TotalBands, NumSubspaces); 
    selectedBands = floor(selectedBands); % Coarse band grouping

    R(1:TotalBands) = inf; % Record the value of each band to the center
    L(1:TotalBands) = 0;

    for iteration = 1:5
        for i = 2:length(selectedBands)
            next = selectedBands;
            class = i - 1;
            index1 = floor((selectedBands(i) - selectedBands(i-1))/2); % The center of each subinterval is obtained as the cluster center
            t = selectedBands(i-1) + index1;  % Current center
            ZL = [];
            ZR = [];
            indexL = [];
            indexR = [];

            if i == 2
               indexr = floor((selectedBands(i+1) - selectedBands(i))/2);
               tr = selectedBands(i) + indexr; % Get the distance to the next center (right)

               ZL = SimilarityMatrix(t, 1:(t-1)); % The left of the class center     
               ZR = SimilarityMatrix(t+1:tr); % The right of the class center   

               indexL = find(ZL < R(1:(t-1)));
               R(indexL) = ZL(indexL);
               L(indexL) = class;

               indexR = find(ZR < R(t+1:tr));
               R(t+indexR) = ZR(indexR);
               L(t+indexR) = class;

            elseif i == length(selectedBands)     
               indexl = floor((selectedBands(i) - selectedBands(i-1))/2);
               tl = selectedBands(i) - indexl; % Get the distance to the next center (left)

               ZL = SimilarityMatrix(t, tl:(t-1)); % The left of the class center       
               ZR = SimilarityMatrix(t, (t+1):TotalBands); % The right of the class center 

               indexL = find(ZL < R(tl:(t-1)));
               R(tl+indexL) = ZL(indexL);       
               L(tl+indexL) = class;

               indexR = find(ZR < R(t+1:TotalBands));
               R(t+indexR) = ZR(indexR);
               L(t+indexR) = class;    

            else
               indexr = floor((selectedBands(i+1) - selectedBands(i))/2);
               tr = selectedBands(i) + indexr; % Get the distance to the last center (right)

               indexl = floor((selectedBands(i-1) - selectedBands(i-2))/2);
               tl = selectedBands(i-1) - indexl; % Get the distance to the last center (left)

               if(tl < 1)
                   tl = 1;
               end

               ZL = SimilarityMatrix(t, tl+1:(t-1)); % The left of the class center 
               ZR = SimilarityMatrix(t+1:tr-1); % The right of the class center     

               indexL = find(ZL < R(tl+1:(t-1)));
               R(tl+1+indexL) = ZL(indexL);       
               L(tl+1+indexL) = class;

               indexR = find(ZR < R(t+1:tr-1));
               R(t+indexR) = ZR(indexR);
               L(t+indexR) = class;
            end   
        end

        [~, t1] = find(L == 0);
        L(t1) = 1;

        %% Remove noise
        for w = 1:NumSubspaces-1
            [~, m] = find(L == w);
            if isempty(m)
               m = selectedBands(w) + 1;
            end

            if m(1) ~= selectedBands(w) + 1
                m = [selectedBands(w) + 1:m(1)-1, m];
            end
            % Compute model
            index = 0;
            for j = 1:length(m) - 1
                if m(j+1) - m(j) > length(m(j+1:end))
                    index = j;
                    break;
                end
            end
            if index == 0
                L(m(1):m(end)) = w;
                selectedBands(w+1) = m(end);
            else
                L(m(index+1:end)) = L(m(index+1)-1);
                L(m(1):m(index)) = w;
                selectedBands(w+1) = m(index);
            end
        end
         if next == selectedBands
             break;
        end

    end
end
