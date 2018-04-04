function [gbih,gbiv,gs6,ga2i] = meteo_sc(rroi,rTi,rpci,freq,rWi,rgai,rmeteo,gbih,gbiv,gs6,ga2i)
%the scattering coefficient of only partly recrystalized snow, linear combination of iborn, wahl==6
[dumgbih,dumgbiv,ags6,dumga2i] = sccoeff(rroi,rTi,rpci,freq,rWi,rgai,6);
gs6=gs6-gs6.*rmeteo+ags6.*rmeteo;

end %endfunction
