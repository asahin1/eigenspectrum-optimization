function [value,gradient] =multiDOFObj(w,params)

    %% Collect parameters
    gridSize = params.gridSize;
    a = params.a;

    %% Collect Matrices and Eigendecompose
    [A,D,L] = generateMatricesFromWeights(w,gridSize);
    L_eps = L + params.eps*eye(length(L));
    [U,diag_lambda] = eig(L_eps);
    lambda = diag(diag_lambda);

    h = a; % without replacement
    n = length(L_eps);
    trL = trace(L_eps);
    trLsqrt = trace(L_eps^(1/2));
    trL2 = trace(L_eps^2);
    b1 = trL/n;
    b2 = sqrt(n*trL2 - trL^2)/n;
    b3 = trLsqrt/n;
    b4 = sqrt(n*trL + n*h^2 - trLsqrt^2)/n;

    %% Compute objective value
    value = h*(b2+2*b3*b4)/(b2*b4*(b2^2+(b1-b3^2)^2+4*b2*b3*b4+2*(b1+b3^2)^2*b4^2+b4^4));

    %% Compute the gradient
    if nargout > 1
        n = length(lambda);
        sumOverL = zeros(length(lambda),1);
        for i=1:length(lambda)
            lambda_i = lambda(i);
            for l=1:length(lambda)
                lambda_l = lambda(l);
                % without h replacement
                if i == l
                    sumOverL(i) = sumOverL(i) + lambda_i.^(-2).*h.*((lambda_i+(-1).*omega_0.^2).^2+2.*(lambda_i+ ...
                                  omega_0.^2).*h.^2+h.^4).^(-2).*((-1).*omega_0.^6+omega_0.^4.*(4.* ...
                                  lambda_i+(-3).*h.^2)+(-1).*(lambda_i+h.^2).^2.*(2.*lambda_i+h.^2)+ ...
                                  (-1).*omega_0.^2.*(lambda_i.^2+3.*h.^4));
                end
                    sumOverL(i) = sumOverL(i) + (-1).*lambda_i.^(-1/2).*lambda_l.^(-1).*omega_0.*h.*((lambda_l+( ...
                                  -1).*omega_0.^2).^2+2.*(lambda_l+omega_0.^2).*h.^2+h.^4).^(-2).*(( ...
                                  -3).*lambda_l.^2+omega_0.^4+(-2).*lambda_l.*h.^2+h.^4+2.* ...
                                  omega_0.^2.*(lambda_l+h.^2)).*n.^(-1);
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

