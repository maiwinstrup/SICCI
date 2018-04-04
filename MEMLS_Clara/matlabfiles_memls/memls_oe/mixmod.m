function [epsi,epsii] = mixmod(f,Ti,Wi,epsi,epsii)

%   calculates the permittivity for Wetness > 0
%      Physical Mixing Model Weise 97 after Matzler 1987 (corrected)
%      water temperature is assumed constant at 273.15 K
%
%   [epsi,epsii] = mixmod(f,Ti,Wi,epsi,epsii)
%       epsi:  real part of the permittivity
%       epsii: imaginary part of the permittivity
%       f:     frequency [GHz]
%       Ti:    physical snow temperature
%       Wi:    wetness [%], no! in vol. frac. 0-1. according to mail 06/01/05 Mätzler -rtt
%       epsi:  real part of dry snow perm.
%       epsii: imaginary part of dry snow perm.
%
%   Version history:
%      1.0    wi 15.7.95
%   
%   Uses: -
%
%
%   Copyright (c) 1997 by the Institute of Applied Physics,
%   University of Bern, Switzerland

Aa = 0.005;
Ab = 0.4975;
Ac = 0.4975;
euw = 4.9;
esw = 88.045; 
frw = 0.11109; % inverse relaxation frequency of water

esa = (esw - epsi)./(3.*(1+Aa.*(esw./epsi-1)));
esb = (esw - epsi)./(3.*(1+Ab.*(esw./epsi-1)));
esc = (esw - epsi)./(3.*(1+Ac.*(esw./epsi-1)));
eua = (euw - epsi)./(3.*(1+Aa.*(euw./epsi-1)));
eub = (euw - epsi)./(3.*(1+Ab.*(euw./epsi-1)));
euc = (euw - epsi)./(3.*(1+Ac.*(euw./epsi-1)));

fa = 1 + Aa * (esw-euw)./(epsi+Aa.*(euw-epsi));
fb = 1 + Ab * (esw-euw)./(epsi+Ab.*(euw-epsi));
fc = 1 + Ac * (esw-euw)./(epsi+Ac.*(euw-epsi));

eea = esa - eua;
eeb = esb - eub;
eec = esc - euc;

fwa = frw./fa;
fwb = frw./fb;
fwc = frw./fc;

depsia = eua + eea ./ (1+(fwa.*f).^2);
depsib = eub + eeb ./ (1+(fwb.*f).^2);
depsic = euc + eec ./ (1+(fwc.*f).^2);
depsi = Wi .* (depsia+depsib+depsic);

depsiia = fwa.*f.*eea ./ (1+(fwa.*f).^2);
depsiib = fwb.*f.*eeb ./ (1+(fwb.*f).^2);
depsiic = fwc.*f.*eec ./ (1+(fwc.*f).^2);
depsii = Wi .* (depsiia + depsiib + depsiic);

epsi = epsi + depsi;
epsii = epsii + depsii;
