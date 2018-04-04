function [gbih,gbiv,gs6,ga2i] = scice_my(si,gbih,gbiv,gs6,ga2i,Ti,dens,freq,pci,sal)

%calculates the scattering coefficient of MY ice from structural parameters

k=(2 .*3.14159) ./(0.3 ./freq);
eice=3.15+0.002*i;
%T=Ti-273.15;
epsi=real(eice);
epsii=imag(eice);

%permittivity of saline ice
[sepsi, sepsii] = sie(si,sal,Ti,freq,epsi,epsii);

eice=sepsi+sepsii*i;

vola=(0.926-dens)./0.926;
emis=eice_s2p(eice,1.0+0.0*i,vola);
ags6=iborn_s2p(eice,1.0+0.0*i,emis,vola,k,pci.*0.001);
gs6=gs6-gs6.*si+ags6.*si;

end %endfunction