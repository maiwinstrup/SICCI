function [pcc]=pccdistr(hs,nsnow,nice)
%compute the scatter size
snow_top=0.07; %the correlation length of new snow [mm]
snow_bottom=0.3; %the correlation length of old snow [mm]
icepcc=1.5; %the correlation length of multiyear ice [mm]

%denne beregning aendrer spredningen for tynd sne, den burde foregå et andet sted

ds=(snow_bottom-snow_top)./nsnow;

nelements=nsnow+nice;
pcc=ones(1,nelements);
ii=10.0/nsnow;
pcc(1)=snow_top;
   for i=2:nsnow
      %pcc(i)=snow_top+(ds.*(i-1));
      %pcc(i)=snow_top+0.0025.*((i).*ii).^2;
      pcc(i)=snow_top+0.001.*((i).*ii).^2;
   end %endfor

   for i=1:nice
      pcc(i+nsnow)=icepcc;
   end %endfor
factor=hs/0.3;
pcc(2:nsnow)=0.05+factor.*(pcc(2:nsnow)-0.05);
%plot(pcc)
end
