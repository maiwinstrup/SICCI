function [gbih,gbiv,gs6,ga2i] = recry_sc(rroi,rTi,rpci,freq,rWi,rgai,rrecry,gbih,gbiv,gs6,ga2i)
%the scattering coefficient of fully recrystalized snow, iborn, wahl==5
[dumgbih,dumgbiv,ags6,dumga2i] = sccoeff(rroi,rTi,rpci,freq,rWi,rgai,5);
gs6=gs6-gs6.*rrecry+ags6.*rrecry;

end %endfunction
