function [ddens]=densdistr(S,Ta,sd,hs,hi,nsnow,nice)
%compute the snow and ice density

[Tsi,dt]=tdistr(Ta,hs,hi,nsnow,nice);
[dS]=Sdistr(S,nsnow,nice);

nelements=nsnow+nice;
ddens=ones(1,nelements);

   for i=1:nsnow
      ddens(i)=sd;
   end %endfor

   for i=1:nice
      ddens(i+nsnow)=dens(dt(i+nsnow),dS(i+nsnow),0.0);
   end %endfor
%plot(ddens)
end
