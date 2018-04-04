function [gbih,gbiv,gs6,ga2i] = scice(si,gbih,gbiv,gs6,ga2i,Ti,sal,freq,pci)

%   calculates the scattering coefficient from structural parameters

k=(2 .*3.14159) ./(0.3 ./freq);
eice=3.15+0.002*i;
T=Ti-273.15;

volb=Vb(T,sal);
salb=Sb(T);
N=Nsw(salb);
conb=condbrine(T,N);
relax=relaxt(T,N);
eb0=epsib0(T,N);
[eb,ebi]=ebrine(freq,eb0,relax,conb); ebri=eb-ebi*i;
%emis=eice_rn2p(eice,eb+ebi*i,volb);
emis=eice_s2p(eice,eb+ebi*i,volb);
ags6=iborn_s2p(eice,ebri,emis,volb,k,pci.*0.001);
gs6=gs6-gs6.*si+ags6.*si;

end %endfunction