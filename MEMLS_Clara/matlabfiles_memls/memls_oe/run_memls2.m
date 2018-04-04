function [amsr]=run_memls2(Ta,S,sd,hs,hi,nsnow,nice)

%Ta=250.0;
%S=0.0;
%hs=0.1;
%hi=1.0;
%nsnow=10;
%nice=20;

[dv,typev,siv,zerov,numv]=ddistr(hs,hi,nsnow,nice);
[Tsiv,Tv]=tdistr(Ta,hs,hi,nsnow,nice);
[Densv]=densdistr(S,Ta,sd,hs,hi,nsnow,nice);
[pccv]=pccdistr2(hs,nsnow,nice);
[dsv]=Sdistr(S,nsnow,nice);

numv=masucolumn(numv);
Tv=masucolumn(Tv);
typev=masucolumn(typev);
Densv=masucolumn(Densv);
dv=masucolumn(dv);
pccv=masucolumn(pccv);
dsv=masucolumn(dsv);
zerov=masucolumn(zerov);
zerov=masucolumn(zerov);
siv=masucolumn(siv);

num=numv;
Ti=flipud(Tv);
typi=flipud(typev);
Densi=flipud(Densv);
di=flipud(dv);
pcci=flipud(pccv);
sali=flipud(dsv);
wci=flipud(zerov);
rmsi=flipud(zerov);
sii=flipud(siv);

[data]=memls_mod(num,di,Ti,wci,Densi,pcci,sali,typi,sii);

%amsr data 10 channels 6v 6h 10v 10h 18v 18h 37v 37h 89v 89h
amsr=data([3:8 11:12 15:16]);

%plot(amsr)
end
