omega = 0:0.001:10;

n = length(L_init);
L_init_eps = L_init+params.eps*eye(n);
sqrtL = (L_init_eps)^(1/2);

pSim = zeros(length(omega),1);
meanSqrtEig = mean(sqrt(lambda_init));
h = params.a;
for i=1:length(pSim)
    pSim(i) = h/(pi*((omega(i)-meanSqrtEig)^2+h^2));
end

h_bar = sqrt(n*h^2 - trace(L_init_eps) + trace(sqrtL)^2/n);

scale = h*n/h_bar;

pObj = zeros(length(omega),1);
for i=1:length(pObj)
    pObj(i) = h_bar*scale/(pi*(norm(omega(i)*eye(n) - sqrtL,"fro")^2+h_bar^2));
end


plot(omega,pSim,"green")
hold on
plot(omega,pObj,"red")
