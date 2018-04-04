function [gai] = abscoeff(epsi,epsii,Ti,freq,Wi)

%   computes the absorption coefficient from the dielectric properties
%
%
%   [gai] = abscoeff(epsi,epsii,Ti,freq,Wi)
%       gai:   absorption coefficient [m^-1]
%       epsi:  real part diel
%       epsii: imaginary part diel
%       Ti:    physical temperature
%       freq:  frequency [GHz]
%       Wi:    volumetric liquid water content
%
%   Version history:
%      1.0    wi 15.7.95
%      1.1    wi 12.11.97 more precise formula for gai used
%   
%   Uses:
%
%
%
%   Copyright (c) 1997 by the Institute of Applied Physics, 
%   University of Bern, Switzerland

% constants
c = 2.99793;
%ny formel for Ka: Ulaby et al. 1981 vol. 1 formel4.140b
lamd=c./(10.*freq);
gai=(4.0.*pi./lamd).*imag(sqrt(epsi+epsii*i));

% Absorptioncoefficient may be suitable for snow but not for saline ice.
%gai = ((2*pi*10*freq).*epsii)./(c.*sqrt(epsi - (epsii.^2./4.*epsi)));
