function entrainment

close all;

num_bins = 50;

%% measure RA entrainment

freqRA = 40;
ampsRA = fliplr([10 12 17 19 99 148]/1000);

a = affpop_single_models();
a.afferents = a.afferents(a.iRA);
edges_dur = linspace(0,0.1,num_bins);
edges_cyc = linspace(0,0.025,num_bins);

figure(1)
set(1,'position',[680   155   654   823], 'PaperPosition', [.25 .25 7.5 7.5]);

c = 1;
for i=1:length(ampsRA)
    s(i) = stim_sine(freqRA,ampsRA(i),[],[],[],20000,[],1);
    r(i) = response(a,s(i));
    
    for j=1:a.num
        if length(r(i).responses(j).spikes)<=1
            continue;
        end
        
        idurRA(j,:) = histc(diff(r(i).responses(j).spikes),edges_dur);
        icycRA(j,:) = histc(mod(r(i).responses(j).spikes,0.025),edges_cyc);
    end
    
    subplot(length(ampsRA),2,c);
    plot(edges_dur,idurRA,'LineWidth',2)
    box off
    c = c+1;
    
    subplot(length(ampsRA),2,c);
    plot(edges_cyc,icycRA,'LineWidth',2)
    box off
    title([num2str(ampsRA(i)*1000) '\mum'])
    c = c+1;
end

screenshot(1,'figs/entrainment01_curr','pdf')
close(1)


%% measure PC entrainment

freqPC = 300;
ampsPC = fliplr([1 5 15 80 100]/1000);

a = affpop_single_models();
a.afferents = a.afferents(a.iPC);
edges_durPC = linspace(0,0.015,num_bins);
edges_cycPC = linspace(0,1/300,num_bins);

figure(2)
set(2,'position',[680   155   654   823], 'PaperPosition', [.25 .25 7.5 7.5]);

c = 1;
for i=1:length(ampsPC)
    s(i) = stim_sine(freqPC,ampsPC(i),[],[],[],20000,[],1);
    r(i) = response(a,s(i));
    
    for j=1:a.num
        idurPC(j,:) = histc(diff(r(i).responses(j).spikes),edges_durPC);
        icycPC(j,:) = histc(mod(r(i).responses(j).spikes,1/300),edges_cycPC);
    end
    
    subplot(length(ampsPC),2,c);
    plot(edges_durPC,idurPC,'LineWidth',2)
    box off
    c = c+1;
    
    subplot(length(ampsPC),2,c);
    plot(edges_cycPC,icycPC,'LineWidth',2)
    box off
    title([num2str(ampsPC(i)*1000) '\mum'])
    c = c+1;
end

screenshot(2,'figs/entrainment02_curr','pdf')
close(2)
