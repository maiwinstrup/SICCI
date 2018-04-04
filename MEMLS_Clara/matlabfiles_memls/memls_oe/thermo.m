function hs=thermo(hi,ks,ki)

%ks=0.30;
%ki=2.1;
hs=ones(1,5)
Ta=[-5,-12.5,-20,-25,-30];
Ta=273.15+Ta;
Tsi=[-3.4,-7.85,-12.3,-15.3,-18.3];
Tsi=273.15+Tsi;
Tw=271.35;
%hi=[0.3,0.5,1.0,2.0]
for i=1:5
%for j=1:4
hs(i)=(hi.*(ks.*Tsi(i)-ks.*Ta(i)))./(ki.*Tw-ki.*Tsi(i))
%end
end
f=(ks.*hi)./(ki.*hs)
Tsi=(Tw+f.*Ta)./(f+1.0)

plot(Ta,hs)

end %endfunction
