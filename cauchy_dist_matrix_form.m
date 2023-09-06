
omega = linspace(0.1,sqrt(lambda_init(end)),10);
pdf_val_exact_init = zeros(length(omega),1);
pdf_val_exact_final = zeros(length(omega),1);
for i=1:length(omega)
    sum_init = 0;
    sum_final = 0;
    for j=1:length(lambda_init)
    % for j=50:50
        sum_init = sum_init + (params.a*sqrt(lambda_init(j))/(pi*((omega(i)-sqrt(lambda_init(j)))^2 + params.a^2*lambda_init(j))));
        sum_final = sum_final + (params.a*sqrt(lambda_final(j))/(pi*((omega(i)-sqrt(lambda_final(j)))^2 + params.a^2*lambda_final(j))));
    end
    pdf_val_exact_init(i) = sum_init;
    pdf_val_exact_final(i) = sum_final;
end

pdf_val_matrix_init = zeros(length(omega),1);
L_init_eps = L_init + params.eps*eye(length(L_init));
pdf_val_matrix_final = zeros(length(omega),1);
L_final_eps = L_final + params.eps*eye(length(L_final));
for i=1:length(omega)
    pdf_val_matrix_init(i) = 1000/((norm(omega(i)*L_init_eps^(-1/2)-eye(length(L_init_eps)),"fro")^2+1));
    pdf_val_matrix_final(i) = 1000/((norm(omega(i)*L_final_eps^(-1/2)-eye(length(L_final_eps)),"fro")^2+1));
end

figure()
subplot(2,1,1)
hist_init = histogram(sqrt(lambda_init),10);
hold on
plot(omega,pdf_val_exact_init,"green")
hold on
plot(omega,pdf_val_matrix_init,"red")
title("Initial graph")
xlabel("Omega")
ylabel("Frequency/Probability")
legend("Square root of eigenvalues","Exact pdf","Approximated pdf")

subplot(2,1,2)
hist_final = histogram(sqrt(lambda_final),10);
hold on
plot(omega,pdf_val_exact_final,"green")
hold on
plot(omega,pdf_val_matrix_final,"red")
title("Final graph")
xlabel("Omega")
ylabel("Frequency/Probability")
legend("Square root of eigenvalues","Exact pdf","Approximated pdf")