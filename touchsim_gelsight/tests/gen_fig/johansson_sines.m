function johansson_sines()

frq=[0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 400];
amp=.5./2.^(0:6); % 0:-6:-48 dB.
pre_indent=.6;
sf=20000;
rad=6;

dur=5./frq;% 5 cycles

ap=affpop_single_models;
nspikes=zeros(length(frq),length(amp),length(ap.afferents));
for ii=1:length(frq)
    for jj=1:length(amp)
        t=(1/sf:1/sf:dur(ii))';
        s=Stimulus(pre_indent+amp(jj)+amp(jj)*sin(2*pi*frq(ii)*t-pi/2),[0 0],sf,rad);
        r=ap.response(s);
        nspikes(ii,jj,:)=dur(ii)*cat(1,r.responses(:).rate);
    end
end
spkpercycle=nspikes/5;

%%
figure('pos',[680   60   400   920])

spkpercycleSA=median(spkpercycle(:,:,ap.iSA1),3);
spkpercycleRA=median(spkpercycle(:,:,ap.iRA),3);
spkpercyclePC=median(spkpercycle(:,:,ap.iPC),3);

ax=[subplot(311) subplot(312) subplot(313)];
h(:,1)=plot(ax(1),frq,spkpercycleSA);
h(:,2)=plot(ax(2),frq,spkpercycleRA);
h(:,3)=plot(ax(3),frq,spkpercyclePC);

set(ax,'xtick',frq,'xscale','log','box','off','linew',1.5,...
    'tickdir','out','xlim',[.4 500],'nextplot','add','colororder',parula)
set(h(:),'linew',2)

for ii=1:3
    xlabel(ax(ii),'Frequency [Hz]')
    ylabel(ax(ii),'Spikes per cycle')
end

title(ax(1),['SA1 (' num2str(length(find(ap.iSA1))) ' units)'])
title(ax(2),['RA (' num2str(length(find(ap.iRA))) ' units)'])
title(ax(3),['PC (' num2str(length(find(ap.iPC))) ' units)'])

screenshot('figs/johansson_sines01_curr','pdf')
