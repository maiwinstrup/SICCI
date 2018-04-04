function [dS]=Sdistr(S,nsnow,nice)
%compute the salinity distribution
%ice_S_bottom=5.0; %the salinity of the bottom ice [ppt]
ice_S_bottom=3.0; %the salinity of the bottom ice [ppt]

dsz=(ice_S_bottom-S)./nice;

nelements=nsnow+nice;

dS=ones(1,nelements);

   for i=1:nsnow
      dS(i)=0.0;
   end %endfor

   for i=1:nice
      dS(i+nsnow)=S+((i-1).*dsz);
   end %endfor
%plot(dS)
end
