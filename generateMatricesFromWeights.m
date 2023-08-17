function [A,D,L] = generateMatricesFromWeights(weightVector,gridSize)

    if(length(weightVector) ~= (gridSize*(gridSize-1)*2))
        error("Weight vector and grid size mismatch!");
    end

    %% Generate the adjacency matrix assuming a uniform grid
    A = zeros(gridSize*gridSize);
    idx = 1;
    for i=1:gridSize
        for j=1:gridSize
            node_idx = (i-1)*gridSize + j;
            if(j ~= gridSize)
                A(node_idx,node_idx+1) = weightVector(idx);
                idx = idx+1;
            end
            if(i ~= gridSize)
                A(node_idx,node_idx+gridSize) = weightVector(idx);
                idx = idx+1;
            end
        end
    end

    % Symmetric elements
    for i=1:length(A)
        for j=i+1:length(A)
            A(j,i) = A(i,j);
        end
    end

    %% Generate Degree and Laplacian Matrices
    D = zeros(size(A));
    for i=1:length(D(:,1))
        D(i,i) = sum(A(i,:));
    end
    L = D - A;
    
end

