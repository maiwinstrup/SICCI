function forward(snow_depth, ice_thickness, temperature, salinity)

ifile='mysnow.dat'
eval (['load ',ifile])
eval (['di = ',ifile(1:6),'(:,5);']);
ii = find(di > 0);
eval (['num = ',ifile(1:6),'(ii,1);']);
eval (['T = ',ifile(1:6),'(ii,2);']);
eval (['typ = ',ifile(1:6),'(ii,3);']);
eval (['Dens = ',ifile(1:6),'(ii,4);']);
eval (['d = ',ifile(1:6),'(ii,5);']);
eval (['pcc = ',ifile(1:6),'(ii,6);']);
eval (['sal = ',ifile(1:6),'(ii,7);']);
eval (['wc = ',ifile(1:6),'(ii,8);']);
eval (['rms = ',ifile(1:6),'(ii,9);']);
eval (['si = ',ifile(1:6),'(ii,10);']);

num=masucolumn(num);
T=masucolumn(T);
typ=masucolumn(typ);
Dens=masucolumn(Dens);
d=masucolumn(d);
pcc=masucolumn(pcc);
sal=masucolumn(sal);
wc=masucolumn(wc);
rms=masucolumn(rms);
si=masucolumn(si);

num=num;
Ti=flipud(T);
typi=flipud(typ);
Densi=flipud(Dens);
di=flipud(d);
pcci=flipud(pcc);
sali=flipud(sal);
wci=flipud(wc);
rmsi=flipud(rms);
sii=flipud(si);

[data]=memls_mod(num,di,Ti,wci,Densi,pcci,sali,typi,sii);

amsr=data([3:12 15:16])

end
