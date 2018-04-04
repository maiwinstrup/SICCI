function [Tsi,T]=tdistr(Ta,hs,hi,nsnow,nice)
%compute the snow ice interface temperature
ki=2.0; %the thermal conductivity of ice [W/mK]
ks=0.31; %the thermal conductivity of snow [W/mK]
%the thermal conductivity is VERY important for the retrieval
Tw=271.35; %the temperature of water under ice

f=(ks.*hi)./(ki.*hs);
Tsi=(Tw+f.*Ta)./(f+1.0);

dTs=(Tsi-Ta)./nsnow;
dTi=(Tw-Tsi)./nice;

nelements=nsnow+nice;
T=ones(1,nelements);

   for i=1:nsnow
      T(i)=Ta+(dTs.*(i-1));
   end %endfor

   for i=1:nice
      T(i+nsnow)=Tsi+(dTi.*(i-1));
   end %endfor
%plot(T)
end
