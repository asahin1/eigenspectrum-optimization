function [value,gradient] =altObjectiveFunctionWithGradient(w,params)

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
        % h = a; % without replacement
        h = a*sqrt(lambda_i); % with replacement
        for j=1:length(lambda)
            lambda_j = lambda(j);
            value = value + ((h*(h^2+lambda_j)+h*lambda_i)...
                /(lambda_j*(lambda_i^2+2*lambda_i*(h^2-lambda_j)+(h^2+lambda_j)^2)));
        end
    end

    %% Compute the gradient
    if nargout > 1
        sumOverL = zeros(length(lambda),1);
        for i=1:length(lambda)
            lambda_i = lambda(i);
            for l=1:length(lambda)
                lambda_l = lambda(l);
                % without h replacement
                % sumOverL(i) = sumOverL(i) + (-(h*(h^6*lambda_l+h^4*(lambda_i+lambda_l)*(lambda_i+3*lambda_l)+...
                %               (lambda_i-lambda_l)*(lambda_i+lambda_l)*(lambda_i^2+4*lambda_i*lambda_l-lambda_l^2)+...
                %               h^2*(2*lambda_i^3+3*lambda_i^2*lambda_l+3*lambda_l^3)))/(lambda_i^2*lambda_l*(h^4+(lambda_i-lambda_l)^2+2*h^2*(lambda_i+lambda_l))^2));
                % with h replacement
                sumOverL(i) = sumOverL(i) + ((1/2).*lambda_i.^(-1/2).*lambda_l.^(-1).*a.*(lambda_l.^2+2.* ...
                              lambda_i.*lambda_l.*((-1)+a.^2)+lambda_i.^2.*(1+a.^2).^2).^(-2).*( ...
                              lambda_l.^3+(-1).*lambda_i.^3.*(1+a.^2).^3+lambda_i.*lambda_l.^2.* ...
                              (5+a.^2)+(-1).*lambda_i.^2.*lambda_l.*(5+6.*a.^2+a.^4))+(-1).*lambda_i.^(-2).*lambda_l.^(1/2).*a.*(lambda_i.^2+2.* ...
                              lambda_i.*lambda_l.*((-1)+a.^2)+lambda_l.^2.*(1+a.^2).^2).^(-2).*( ...
                              2.*lambda_i.^3+lambda_l.^3.*(1+a.^2).^3+lambda_i.^2.*(lambda_l+5.* ...
                              lambda_l.*a.^2)+4.*lambda_i.*lambda_l.^2.*((-1)+a.^4)));
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

