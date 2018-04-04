function [epsi, epsii] = mysie(si,rho,Ti,sal,freq,epsi,epsii)

% computes the dielectric constant of ice if it is an ice layer
% Background information: Ulaby et al. 1986 vol. III.
% i.e. if si == 1
% si: sea ice/snow [1/0]
% rho: density of icelayer [kg/m3]
% Ti: Thermometric temperature of layer [K]
% sal: salinity of ice [psu]
% freq: Frequency in GHz
% epsi: initial permittivity (of snow)
% epsii: initial loss (of snow)

%permittivity of saline ice
[sepsi, sepsii] = sie(si,sal,Ti,freq,epsi,epsii);

eice=epice(Ti,freq);                            %fresh ice dielectric const
vola=(0.926 - rho)./0.926;                      %volume of air
emis=eice_s2p(sepsi+i*sepsii,1.0+0.0i,vola);    %dielectric constant of sea ice (spherical inclusions)

  aepsi=real(emis);
  aepsii=imag(emis);

  epsi=epsi-epsi.*si+aepsi.*si;
  epsii=epsii-epsii.*si+aepsii.*si;

end %endfunction
