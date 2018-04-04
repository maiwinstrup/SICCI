function [epsi, epsii] = sie(si,sal,Ti,freq,epsi,epsii)

% computes the dielectric constant of ice if it is an ice layer
% Background information: Ulaby et al. 1986 vol. III.
% i.e. if si == 1
% si: sea ice/snow [1/0]
% sal: salinity [ppt/psu]
% Ti: Thermometric temperature of layer [K]
% freq: Frequency in GHz
% epsi: initial permittivity (of snow)
% epsii: initial loss (of snow)

T=Ti-273.15;                          %get therm. temp. in C
eice=epice(Ti,freq);                  %fresh ice dielectric const
volb=Vb(T,sal);                       %volume of brine
salb=Sb(T);                           %salinity of brine
N=Nsw(salb);                          %normality of brine solution
conb=condbrine(T,N);                  %conductivity of brine solution
relax=relaxt(T,N);                    %relaxation time of brine
eb0=epsib0(T,N);                      %static dielectric const. of brine
[eb,ebi]=ebrine(freq,eb0,relax,conb); %dielectric constant of brine
emis=eice_rn2p(eice,eb+ebi*i,volb);  %dielectric constant of sea ice (random needles)
%emis=eice_s2p(eice,eb+ebi*i,volb);    %dielectric constant of sea ice (spherical inclusions)
%emis=0.5.*eice_rn2p(eice,eb+ebi*i,volb)+0.5.*eice_s2p(eice,eb+ebi*i,volb); %mix of spheres and random needles

  aepsi=real(emis);
  aepsii=imag(emis);

  epsi=epsi-epsi.*si+aepsi.*si;
  epsii=epsii-epsii.*si+aepsii.*si;
