
com=ATTC;
pil=inr;
maxfb=size(com);
minfb=size(pil);
diff=[];
pil=resample(pil,maxfb,minfb);
diff=com(:,4)-pil(:,3);

maxfb=size(outr);
minfb=size(diff);

if maxfb>minfb
    diff=resample(diff,maxfb,minfb);
else
    outr=resample(outr,minfb,maxfb);
end

outr(:,11)=diff;