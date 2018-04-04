function eice = epsice(Ti,freq)

%   calculates the dielectric permittivity of ice
%   After Hufford, Mitzima and Matzler
%
%   eice = epsice(Ti,freq)
%      eice:  dielectric permittivity of ice
%      Ti:    temperature in K 
%      freq:  frequency in GHz
%
%   Version history:
%      1.0    wi 15.7.95
%   
%   Uses:
%
%
%
%   Copyright (c) 1997 by the Institute of Applied Physics,
%   University of Bern, Switzerland

pp = (300 ./ Ti) -1;
B1 = 0.0207; b = 335.25; B2 = 1.16e-11; db = exp(-10.02+0.0364 .* (Ti-273));
beta = ((B1.* exp(b./Ti)) ./ (Ti .* (exp(b./Ti)-1).^2)) + B2 .* freq^2 + db;
alpha = ((0.00504 + 0.0062.*pp) .* exp(-22.1 .* pp));
eice = (alpha./freq) + (beta.*freq);
