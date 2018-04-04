function se=epsib0(T,N)

% static dielectric constant of brine, se. Ulaby et al. 1986, E65a

% T: thermometric temperature [C]
% N: normality of brine

eps=88.045-0.4147 .*T+6.295E-4 .*(T .^2)+1.075E-5 .*(T .^3);
a=1.-0.255 .*N+5.15E-2 .*(N .^2)-6.89E-3 .*(N .^3);
se=eps.*a;
