function adaptation_prop()
close all;

sf = 20000;

load data/kni_SA
spikes_SA=spikes;

load data/kni_RA
spikes_RA=spikes;


%% SA
indentSA=[1.6:.2:3.2;.4:.2:2;.4:.2:2;.5 .9 1.3 1.8 2.1 2.6  2.9 3.2 3.5]';
aSA = affpop_single_models();
aSA.afferents(~aSA.iSA1)=[];

cols=[252,187,161
252,146,114
251,106,74
239,59,44
203,24,29
153,0,13]/255;

figure(1)
set(1,'pos',[100 400 650 500]);
ax=[subplot(3,2,1) subplot(3,2,3) subplot(3,2,5)...
    subplot(3,2,2) subplot(3,2,4) subplot(3,2,6)];
set(ax,'nextplot','add','colororder',cols);

for ii=1:2:length(indentSA)
    s=stim_ramp(indentSA(ii,3),1,[0 0],sf,indentSA(ii,3)/40,'lin',1);
    r = aSA.response(s);r={r.responses.spikes};
    r=cellfun(@(x) x+.025,r,'uni',0);
    axes(ax(1)); plot(1/sf:1/sf:s.duration,s.trace)
    plot_spikes(r(1),'par',ax(3),'neuron_o',(ii+1)/2-1,'col',cols((ii+1)/2,:),'hold','on')
    
    s=stim_ramp(indentSA(ii,2),1,[0 0],sf,indentSA(ii,2)/40,'lin',1);
    r = aSA.response(s);r={r.responses.spikes};
    axes(ax(4)); plot(1/sf:1/sf:s.duration,s.trace)
    r=cellfun(@(x) x+.025,r,'uni',0);
    plot_spikes(r(3),'par',ax(6),'neuron_o',(ii+1)/2-1,'col',cols((ii+1)/2,:),'hold','on')
    
    plot_spikes(spikes_SA(length(indentSA)+1-ii,3),'neuron_o',(ii+1)/2-1,'par',ax(2),'col',cols((ii+1)/2,:),'hold','on')
    plot_spikes(spikes_SA(length(indentSA)+1-ii,2),'neuron_o',(ii+1)/2-1,'par',ax(5),'col',cols((ii+1)/2,:),'hold','on')
end
for ii=1:6,xlabel(ax(ii),''); end
set(ax,'xlim',[-.01 1.01],'xcolor','none')
%scalebar(ax(3),'x','100 ms',.1,0)
%scalebar(ax(6),'x','100 ms',.1,0)

title(ax(1),'Example 1')
title(ax(4),'Example 2')
ylabel(ax(1),'Depth [mm]')
ylabel(ax(4),'Depth [mm]')
ylabel(ax(2),'Spikes')
ylabel(ax(5),'Spikes')
ylabel(ax(3),'Spikes')
ylabel(ax(6),'Spikes')

for ii=1:3
    set(ax(ii),'pos',get(ax(ii),'pos')+[-.06 0 .06 0])
end
for ii=4:6
    set(ax(ii),'pos',get(ax(ii),'pos')+[0.01 0 .06 0])
end

set(ax([2 3 5 6]),'ycolor','none')

screenshot(1,'figs/adaptation_prop01_curr','pdf')
screenshot(1,'figs/adaptation_prop01_curr','png')
%close(1)

%% RA
indentRA=[1 1 1 1 1 1;2 2 2 2 2 2]';
delaysRA=cellfun(@(x) find(x>.031,1),posall);

for ii=1:length(spikes_RA(:))
    spikes_RA{ii}=spikes_RA{ii}-delaysRA(ii)/1000;
end
spikes_RA{1,2}=spikes_RA{1,2}-.008;
aRA = affpop_single_models();
aRA.afferents(~aRA.iRA)=[];

cols=[218,218,235
188,189,220
158,154,200
128,125,186
106,81,163
74,20,134]/255;

figure(2)
set(2,'pos',[1000 400 650 500]);
ax=[subplot(3,2,1) subplot(3,2,3) subplot(3,2,5)...
    subplot(3,2,2) subplot(3,2,4) subplot(3,2,6)];
set(ax,'nextplot','add','colororder',cols);

for ii=1:length(speeds)
    s=stim_ramp(indentRA(ii,1),2,[0 0],sf,indentRA(ii,1)/speeds(ii,1),'lin',1);
    r = aRA.response(s); r={r.responses.spikes};
    r=cellfun(@(x) x+.0,r,'uni',0);
    axes(ax(1)); plot(1/sf:1/sf:s.duration,s.trace)
    plot_spikes(r(2),'par',ax(3),'neuron_o',ii-1,'col',cols(ii,:),'hold','on')
    
    s=stim_ramp(indentRA(ii,2),2,[0 0],sf,indentRA(ii,2)/speeds(ii,2),'lin',1);
    r = aRA.response(s);r={r.responses.spikes};
    r=cellfun(@(x) x+.005,r,'uni',0);
    axes(ax(4)); plot(1/sf:1/sf:s.duration,s.trace)
    plot_spikes(r(8),'par',ax(6),'neuron_o',ii-1,'col',cols(ii,:),'hold','on')
    
    plot_spikes(spikes_RA(ii,1),'par',ax(2),'neuron_o',ii-1,'col',cols(ii,:),'hold','on')
    plot_spikes(spikes_RA(ii,2),'par',ax(5),'neuron_o',ii-1,'col',cols(ii,:),'hold','on')
end

for ii=1:6,xlabel(ax(ii),''); end
set(ax(1:3),'xlim',[-.001 .101],'xcolor','none')
set(ax(4:6),'xlim',[-.001 .401],'xcolor','none')
%scalebar(ax(3),'x','10 ms',.01,0.01)
%scalebar(ax(6),'x','50 ms',.05,0.01)

title(ax(1),'Example 1')
title(ax(4),'Example 2')
ylabel(ax(1),'Depth [mm]')
ylabel(ax(4),'Depth [mm]')
ylabel(ax(2),'Spikes (actual)')
ylabel(ax(5),'Spikes (actual)')
ylabel(ax(3),'Spikes (model)')
ylabel(ax(6),'Spikes (model)')


set(ax([2 3 5 6]),'ycolor','none')
for ii=1:3
    set(ax(ii),'pos',get(ax(ii),'pos')+[-.06 0 .06 0])
end
for ii=4:6
    set(ax(ii),'pos',get(ax(ii),'pos')+[0.01 0 .06 0])
end

screenshot(2,'figs/adaptation_prop02_curr','pdf')
screenshot(2,'figs/adaptation_prop02_curr','png')
%close(2)