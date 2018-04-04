function N = Nsw(Ssw)

% normality of sea water or brine Ulaby et al. 1986, E20
% Ssw: salinity of brine or sea water

N = 0.91141 .*Ssw .*(1.707e-2 +1.205e-5 .*Ssw+4.058e-9 .*(Ssw .^2));
