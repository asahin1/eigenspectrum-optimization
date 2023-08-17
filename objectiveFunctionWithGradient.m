function [value,gradient] = objectiveFunctionWithGradient(w,params)

    %% Collect parameters
    gridSize = params.gridSize;
    a = params.a;

    %% Collect Matrices and Eigendecompose
    [A,D,L] = generateMatricesFromWeights(w,gridSize);
    [U,diag_lambda] = eig(L + params.eps*eye(length(L)));
    lambda = diag(diag_lambda);

    %% Compute objective value
    value = 0;
    for i=1:length(lambda)
        lambda_i = lambda(i);
        h = a*sqrt(lambda_i);
        for j=1:length(lambda)
            lambda_j = lambda(j);
            denominator = (sqrt(lambda_i)-sqrt(lambda_j))^2 + h*h;
            value = value + (4*h*h/lambda_j)*(1/(denominator*denominator));
        end
    end

    %% Compute the gradient
    if nargout > 1
        sumOverL = zeros(length(lambda),1);
        for i=1:length(lambda)
            lambda_i = lambda(i);
            for l=1:length(lambda)
                lambda_l = lambda(l);
                sumOverL(i) = sumOverL(i) + (-4*a*a*(a*a*lambda_i+lambda_i-lambda_l)) / (lambda_l * (a*a*lambda_i+(sqrt(lambda_i)-sqrt(lambda_l))^2)^3)...
                    + (4*a*a*lambda_l*(-a*a*lambda_l-lambda_l+4*sqrt(lambda_l)*sqrt(lambda_i)-3*lambda_i)) / (lambda_i*lambda_i * (a*a*lambda_l+(sqrt(lambda_l)-sqrt(lambda_i))^2)^3);
            end
        end
        idx = 1;
        gradient = zeros(length(w),1);
        for i=1:gridSize
            for j=1:gridSize
                node_idx = (i-1)*gridSize + j;
                if(j ~= gridSize)
                    % Compute gradient wrt w_(node_idx, node_idx+1)
                    for k=1:length(lambda)
                        gradient(idx) = gradient(idx) + sumOverL(k)*(U(node_idx,k)-U(node_idx+1,k))^2;
                    end
                    idx = idx+1;
                end
                if(i ~= gridSize)
                    % Compute gradient wrt w_(node_idx,node_idx+gridSize)
                    for k=1:length(lambda)
                        gradient(idx) = gradient(idx) + sumOverL(k)*(U(node_idx,k)-U(node_idx+gridSize,k))^2;
                    end
                    idx = idx+1;
                end
            end
        end
    end
end

