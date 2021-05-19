% test_streaming.m

stim = [0 0.25 0.5 0.75 1 1 1 1 1 1 1 1 0.75 0.5 0.25 0];

aSA = AfferentStream('SA1','idx',1);
aRA = AfferentStream('RA','idx',1);
aPC = AfferentStream('PC','idx',1);

for s=1:length(stim)
    rSA(s) = aSA.response(stim(s),40);
    rRA(s) = aRA.response(stim(s),40);
    rPC(s) = aPC.response(stim(s),40);
end

%%

figure
hold on
plot(stim*50,'k')
plot(rSA,'g')
plot(rRA,'b')
plot(rPC,'r')
