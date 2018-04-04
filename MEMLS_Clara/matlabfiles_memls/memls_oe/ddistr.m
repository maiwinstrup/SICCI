function [d,type,si,zero,num]=ddistr(hs,hi,nsnow,nice)
%compute the snow and ice layer thickness d

ds=(hs)./nsnow;
di=(hi)./nice;

nelements=nsnow+nice;
d=ones(1,nelements);
type=ones(1,nelements);
si=ones(1,nelements);
zero=zeros(1,nelements);
num=1:nelements;

   for i=1:nsnow
      d(i)=ds;
      type(i)=2;
      si(i)=0;
   end %endfor

   for i=1:nice
      d(i+nsnow)=di;
      type(i+nsnow)=4;
      si(i+nsnow)=1;
   end %endfor
%plot(d)
end
