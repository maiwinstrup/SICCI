function sigma = condbrine(T,N)

% conductivity of brine Ulaby et al. 1986, vol. III E68,E69
% T: thermometric temperature [C]
% N: normality of brine

D=25. -T;
sig=N .*(10.39-2.378 .*N+0.683 .*(N .^2)-0.135 .*(N .^3)+(1.01e-2) .*(N .^4));
c=1.-(1.96e-2) .*D+(8.08e-5) .*(D .^2)-N .*D .*(3.02e-5)+(3.92e-5) .*D+N .*((1.72e-5)-(6.58e-6) .*D);
sigma=c .*sig;
