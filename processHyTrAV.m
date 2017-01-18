function[] = processHyTrAV()
load('C:\Users\Navaneeth Krishnan B\Desktop\IISc\DataAquisition\11_05_09.bin-1077528.mat');
timeus = (TIME(:,2));
timeus = timeus - (min(TIME(:,2)));

maxsize = size((ATT(:,7)),1);
minsize = size((RC(:,3)),1);
maxsizealt = size(SENS(:,2),1);

gndalt = min(SENS(:,3));
altituderel = (SENS(:,3) - gndalt);
maxaltitude = max(altituderel);
altitudenorm = altituderel/maxaltitude;

altitude = resample(altitudenorm, maxsize,maxsizealt);
inpitch = resample((RC(:,2)), maxsize,minsize);
inroll = resample((RC(:,3)), maxsize,minsize);
inyaw = resample((RC(:,5)), maxsize,minsize);
inthrottle = resample((RC(:,4)), maxsize,minsize);

outpitch = (ATT(:,7));
outroll = (ATT(:,6));
outyaw = (ATT(:,8));

outpitchrate = (ATT(:,10));
outrollrate = (ATT(:,9));
outyawrate = (ATT(:,11));

inpitchc = (ATTC(:,3));
inrollc = (ATTC(:,2));
inyawc = (ATTC(:,4));

outpitchinv = -1 * outpitch;
outpitchrateinv = -1 * outpitchrate;

plot(timeus,outpitchinv,'DisplayName','outpitch');
hold all;
%plot(outpitchrateinv,'DisplayName','outpitchrate');
%plot(outroll,'DisplayName','outroll');
%plot(outrollrate,'DisplayName','outrollrate');
%plot(altitude,'DisplayName','altitude');
plot(timeus,inpitch,'DisplayName','inpitchc');
plot(timeus,inroll,'DisplayName','inrollc');
plot(timeus,inyaw,'DisplayName','inyawc');
%plot(inthrottle,'DisplayName','throttle');
hold off;

filename = max(timeus);
fname = sprintf('SysID%d.mat', filename);
save(fname)