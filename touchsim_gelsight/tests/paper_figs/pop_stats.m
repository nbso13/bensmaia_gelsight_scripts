% pop_stats.m

a = affpop_hand('D2');
a2 = affpop_hand('D1');
a3 = affpop_hand();
a.afferents = [a.afferents a2.afferents];
a.afferents(a.iPC) = [];
a.afferents = [a.afferents a3.afferents(a3.iPC)];
clear a2;

aSA = AfferentPopulation;
aSA.afferents = a.afferents(a.iSA1);
aRA = AfferentPopulation;
aRA.afferents = a.afferents(a.iRA);
aPC = AfferentPopulation;
aPC.afferents = a.afferents(a.iPC);

fprintf('Generated afferent population.\n')

%% simulate vibrotactile stimuli

freqs = [15  300];
amps =  [500 50];

% calculate responses from peripheral model
tic
for i=1:numel(freqs)
    s = stim_sine(freqs(i),amps(i)/1000,[],[],[],5000,[],0.5);
    r_freq(i) = a.response(s,1);
    rates_freq(:,i) = r_freq(i).rate;
    
    lats_freq(:,i) = NaN*zeros(size(rates_freq,1),1);
    ll_tmp = {r_freq(i).responses(:).spikes};
    lats_freq(rates_freq(:,i)>0,i) = cellfun(@(x) x(1),ll_tmp(rates_freq(:,i)>0));
end
toc
fprintf('Calculated response to skin vibrations.\n')

clear s

ratesSA = rates_freq(a.iSA1,:);
ratesRA = rates_freq(a.iRA,:);
ratesPC = rates_freq(a.iPC,:);

ratesSA(ratesSA==0) = NaN;
ratesRA(ratesRA==0) = NaN;
ratesPC(ratesPC==0) = NaN;

rSA = ratesSA;
rRA = ratesRA;
rPC = ratesPC;

nSA = ~isnan(ratesSA);
nRA = ~isnan(ratesRA);
nPC = ~isnan(ratesPC);

latsSA = lats_freq(a.iSA1,:);
latsRA = lats_freq(a.iRA,:);
latsPC = lats_freq(a.iPC,:);

%% simulate grasp stimuli

tic
contacts = [1 3];
[origin,theta,pxl_per_mm,regionprop] = plot_hand(NaN);
rot = [cos(-theta) -sin(-theta);sin(-theta) cos(-theta)];
for c=1:length(contacts)
    contact_locs(c,:) = regionprop(contacts(c)).Centroid;
    
end
contact_locs = bsxfun(@minus,contact_locs,origin);
contact_locs = contact_locs/pxl_per_mm;
contact_locs = contact_locs*inv(rot);

s = pad(pad(stim_indent_shape(contact_locs,stim_ramp(0.75,.5,[],[],0.05,'lin',5)),[0 0.01]),0.05);
r = a.response(s,1);

rr_tmp = r.psth(125)';

rates(:,1) = rr_tmp(:,1)/0.125;
rates(:,2) = sum(rr_tmp(:,2:3),2)/0.25;
toc
fprintf('Calculated response to grasp stimuli.\n')

ratesSA = rates(a.iSA1,:);
ratesRA = rates(a.iRA,:);
ratesPC = rates(a.iPC,:);

ratesSA(ratesSA==0) = NaN;
ratesRA(ratesRA==0) = NaN;
ratesPC(ratesPC==0) = NaN;

rSA(:,end+1:end+2) = ratesSA;
rRA(:,end+1:end+2) = ratesRA;
rPC(:,end+1:end+2) = ratesPC;

nSA(:,end+1:end+2) = ~isnan(ratesSA);
nRA(:,end+1:end+2) = ~isnan(ratesRA);
nPC(:,end+1:end+2) = ~isnan(ratesPC);

%% cleanup

nSA = double(nSA);
nRA = double(nRA);
nPC = double(nPC);

nSA(nSA==0) = NaN;
nRA(nRA==0) = NaN;
nPC(nPC==0) = NaN;

%% export grasp video

fprintf('Exporting videos...\n')

tic
video(r_freq(2),10);
toc

tic
video(r,10);
toc

%% make some figures

% response rates
figure(1)
set(gcf,'pos',[0 0 900 200]);
for i=1:size(rSA,2)
    subplot(1,size(rSA,2),i)
    semilogy(1,nansum(rSA(:,i)),'o','MarkerSize',7,'LineWidth',4,'color',affcol(1))
    hold on
    semilogy(2,nansum(rRA(:,i)),'o','MarkerSize',7,'LineWidth',4,'color',affcol(2))
    semilogy(3,nansum(rPC(:,i)),'o','MarkerSize',7,'LineWidth',4,'color',affcol(3))
    
    xlim([.5 3.5])
    ylim([1 60000])
    box off
    set(gca,'xtick',[1 2 3],'xticklabel',{'SA1','RA','PC'},'ytick',[1 100 10000])
end
screenshot(1,'pop_stats1','pdf')
close(1)

% spatial SA
figure(2)
set(gcf,'pos',[0 0 900 250]);
for i=1:size(rSA,2)
    subplot(1,size(rSA,2),i)
    col = repmat([.8 .8 .8],aSA.num,1) - repmat(rSA(:,i),1,3).*repmat(([.8 .8 .8]-affcol(1))./max(rSA(:,i)),aSA.num,1);
    plot_hand(gca,'names',false,'axes',false,'centers',false,'scalebar',false,'region','D2',...
        'afferents',aSA,'color',col);
end
screenshot(2,'pop_stats2','pdf')
close(2)

% spatial RA
figure(3)
set(gcf,'pos',[0 0 900 250]);
for i=1:size(rSA,2)
    subplot(1,size(rSA,2),i)
    col = repmat([.8 .8 .8],aRA.num,1) - repmat(rRA(:,i),1,3).*repmat(([.8 .8 .8]-affcol(2))./max(rRA(:,i)),aRA.num,1);
    plot_hand(gca,'names',false,'axes',false,'centers',false,'scalebar',false,'region','D2',...
        'afferents',aRA,'color',col);
end
screenshot(3,'pop_stats3','pdf')
close(3)

% spatial PC
figure(4)
set(gcf,'pos',[0 0 900 250]);
for i=1:size(rSA,2)
    subplot(1,size(rSA,2),i)
    col = repmat([.8 .8 .8],aPC.num,1) - repmat(rPC(:,i),1,3).*repmat(([.8 .8 .8]-affcol(3))./max(rPC(:,i)),aPC.num,1);
    plot_hand(gca,'names',false,'axes',false,'centers',false,'scalebar',false,...
        'afferents',aPC,'color',col);
end
screenshot(4,'pop_stats4','pdf')
close(4)

% spatial SA thumb
figure(5)
set(gcf,'pos',[0 0 900 250]);
for i=1:size(rSA,2)
    subplot(1,size(rSA,2),i)
    col = repmat([.8 .8 .8],aSA.num,1) - repmat(rSA(:,i),1,3).*repmat(([.8 .8 .8]-affcol(1))./max(rSA(:,i)),aSA.num,1);
    plot_hand(gca,'names',false,'axes',false,'centers',false,'scalebar',false,'region','D1',...
        'afferents',aSA,'color',col);
end
screenshot(5,'pop_stats5','pdf')
close(5)

% spatial RA thumb
figure(6)
set(gcf,'pos',[0 0 900 250]);
for i=1:size(rSA,2)
    subplot(1,size(rSA,2),i)
    col = repmat([.8 .8 .8],aRA.num,1) - repmat(rRA(:,i),1,3).*repmat(([.8 .8 .8]-affcol(2))./max(rRA(:,i)),aRA.num,1);
    plot_hand(gca,'names',false,'axes',false,'centers',false,'scalebar',false,'region','D1',...
        'afferents',aRA,'color',col);
end
screenshot(6,'pop_stats6','pdf')
close(6)

%% collated response rates figure

rr = [nansum(rSA);nansum(rRA);nansum(rPC)];
rr(rr==0) = 1;

figure(7)
set(gcf,'pos',[0 0 200 200]);
semilogy([1 2 3],rr,'-k')
hold on
semilogy(1,rr(1,:),'o','MarkerSize',7,'LineWidth',4,'color',affcol(1))
semilogy(2,rr(2,:),'o','MarkerSize',7,'LineWidth',4,'color',affcol(2))
semilogy(3,rr(3,:),'o','MarkerSize',7,'LineWidth',4,'color',affcol(3))

xlim([.5 3.5])
ylim([1 60000])
box off
set(gca,'xtick',[1 2 3],'xticklabel',{'SA1','RA','PC'},'ytick',[1 100 10000])
ylabel('Total spikes / s')

screenshot(7,'pop_stats7','pdf')
close(7)
