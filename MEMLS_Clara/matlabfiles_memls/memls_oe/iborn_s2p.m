% improved born approximation by C. M?tzler (1998). J. Appl. Phys. 83(11),6111-7

function ss=iborn_s2p(e1,e2,eeff,v,k,pcc)
%scattering coefficient of a collection of spherical
%particles with correlation length pcc using improved born approximation
ss=(3 .*pcc.^3 .*k.^4 ./32) .*v .*(1-v) .*abs(((e2-e1) ...
.*(2 .*eeff+e1)) ./(2 .*eeff+e2)).^2;

end %endfunction
