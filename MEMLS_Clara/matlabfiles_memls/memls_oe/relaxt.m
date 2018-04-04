function relax=relaxt(T,N)

% relaxation time of brine, Ulaby et al. 1986, E67
% T: thermometric temperature
% N: normality of brine

rel=(1.1109E-10)-(3.824E-12) .*T+(6.938E-14) .*(T .^2)-(5.096E-16) .*(T .^3);
b=1.-(0.146E-2) .*T .*N-(4.89E-2) .*N-(2.97E-2) .*(N .^2)+(5.64E-3) .*(N .^3);

relax=rel .*b;
