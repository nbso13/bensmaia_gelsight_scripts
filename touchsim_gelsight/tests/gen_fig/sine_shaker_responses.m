function sine_shaker_responses


%% set up stimuli
freqs = [1 5 10 25 50 100 200 300 400 600 800 1000];
freqs = repmat(freqs,14,1);

amps(:,1) = [5.0000    6.9477    9.6541   13.4148   18.6405   25.9017   35.9915   50.0118   69.4934   96.5641  134.1799  186.4486  259.0781  360.0000]';
amps(:,2) = [5.0000    6.9477    9.6541   13.4148   18.6405   25.9017   35.9915   50.0118   69.4934   96.5641  134.1799  186.4486  259.0781  360.0000]';
amps(:,3) = [2.5000    3.3880    4.5914    6.2222    8.4322   11.4272   15.4861   20.9866   28.4408   38.5427   52.2327   70.7853   95.9275  130.0000]';
amps(:,4) = [2.5000    3.3880    4.5914    6.2222    8.4322   11.4272   15.4861   20.9866   28.4408   38.5427   52.2327   70.7853   95.9275  130.0000]';
amps(:,5) = [1.4835    2.2007    3.2647    4.8431    7.1847   10.6583   15.8114   23.4559   34.7963   51.6196   76.5765  113.5997  168.5228  250.0000]';
amps(:,6) = [0.7671    1.1768    1.8053    2.7696    4.2489    6.5184   10.0000   15.3413   23.5355   36.1064   55.3918   84.9781  130.3673  200.0000]';
amps(:,7) = [0.2015    0.3248    0.5236    0.8440    1.3606    2.1933    3.5355    5.6993    9.1872   14.8097   23.8731   38.4833   62.0349  100.0000]';
amps(:,8) = [0.1008    0.1624    0.2618    0.4220    0.6803    1.0966    1.7678    2.8496    4.5936    7.4048   11.9366   19.2417   31.0175   50.0000]';
amps(:,9) = [0.1918    0.2942    0.4513    0.6924    1.0622    1.6296    2.5000    3.8353    5.8839    9.0266   13.8480   21.2445   32.5918   50.0000]';
amps(:,10) = [0.1250    0.1825    0.2665    0.3890    0.5680    0.8293    1.2108    1.7678    2.5810    3.7683    5.5018    8.0327   11.7279   17.1230]';
amps(:,11) = [0.1250    0.1737    0.2413    0.3353    0.4659    0.6474    0.8996    1.2500    1.7369    2.4134    3.3534    4.6595    6.4743    8.9961]';
amps(:,12) = [0.1250    0.1737    0.2413    0.3353    0.4659    0.6474    0.8996    1.2500    1.7369    2.4134    3.3534    4.6595    6.4743    8.9961]';

%% calculate responses from peripheral model

a = affpop_single_models();

for i=1:numel(freqs)
    s(i) = stim_sine(freqs(i),amps(i)/1000,[],[],[],20000,[],0.5);
    r(i) = response(a,s(i));
    rates(:,i) = r(i).rate;
end

%% calculate absolute and tuning thresholds

thres_num_spikes = 10;

thres = nan(a.num,size(freqs,2));
tuning = nan(a.num,size(freqs,2));
for s=1:size(freqs,2)
    ind_freq = find(freqs==freqs(1,s));
    for aff=1:a.num
        idx = find(rates(aff,ind_freq)>=thres_num_spikes,1);
        if ~isempty(idx)
            thres(aff,s) = amps(idx,s);
        end
        idx = find(rates(aff,ind_freq)>0.8*freqs(1,s),1);
        if ~isempty(idx)
            tuning(aff,s) = amps(idx,s);
        end
    end
    thresRA(s) = mean(denan(thres(a.iRA,s)));
    thresPC(s) = mean(denan(thres(a.iPC,s)));
    tuningRA(s) = mean(denan(tuning(a.iRA,s)));
    tuningPC(s) = mean(denan(tuning(a.iPC,s)));
end

%% fit rate-intensity functions to mean responses

rec_log_lin =@(params,x) rectify(params(1)*(log10(x)-params(2)));
options = optimoptions('lsqcurvefit','Display','off');

for f=1:size(freqs,2)
    pSA(f,:) = lsqcurvefit(rec_log_lin,[50 1],amps(:,f),mean(rates(a.iSA1,(f-1)*14+1:(f-1)*14+14),1)',[],[],options);
    pRA(f,:) = lsqcurvefit(rec_log_lin,[50 1],amps(:,f),mean(rates(a.iRA,(f-1)*14+1:(f-1)*14+14),1)',[],[],options);
    pPC(f,:) = lsqcurvefit(rec_log_lin,[50 1],amps(:,f),mean(rates(a.iPC,(f-1)*14+1:(f-1)*14+14),1)',[],[],options);
end

%% plot average rate-intensity curves (see Muniak et al., Fig. 4)

close all
figure(1)
set(gcf,'pos',[680    49   480   948], 'PaperPosition', [.25 .25 5 10]);
subplot(311)
for f=1:size(freqs,2)
    semilogx(amps(:,f),mean(rates(a.iSA1,(f-1)*14+1:(f-1)*14+14),1),'LineWidth',2)
    hold all
end
title('SA1')
legend(num2str(freqs(1,1:7)'),'Location','NorthWest')

subplot(312)
for f=1:size(freqs,2)
    semilogx(amps(:,f),mean(rates(a.iRA,(f-1)*14+1:(f-1)*14+14),1),'LineWidth',2)
    hold all
end
title('RA')

subplot(313)
for f=1:size(freqs,2)
    semilogx(amps(:,f),mean(rates(a.iPC,(f-1)*14+1:(f-1)*14+14),1),'LineWidth',2)
    hold all
end
title('PC')

screenshot(1,'figs/sine_shaker_responses01_curr','pdf')
close(1)

%% compare parameters of mean rate-intensity functions

figure(2)
set(gcf,'pos',[680    49   240   474]);

subplot(211)
loglog(freqs(1,:),pSA(:,1),'Color',[96 173 172]/255,'LineWidth',2);
hold on
loglog(freqs(1,:),pRA(:,1),'Color',[131 50 132]/255,'LineWidth',2);
loglog(freqs(1,:),pPC(:,1),'Color',[162 121 70]/255,'LineWidth',2);
ylim([1e0 1e3])
xlim([1 1000])
ylabel('\alpha')
box off

subplot(212)
loglog(freqs(1,:),10.^pSA(:,2),'Color',[96 173 172]/255,'LineWidth',2);
hold on
loglog(freqs(1,:),10.^pRA(:,2),'Color',[131 50 132]/255,'LineWidth',2);
loglog(freqs(1,:),10.^pPC(:,2),'Color',[162 121 70]/255,'LineWidth',2);
ylim([1e-1 1e2])
xlim([1 1000])
legend('SA','RA','PC','Location','SouthWest')
xlabel('Frequency (Hz)')
ylabel('10^\beta (\mum)')
box off

screenshot(2,'figs/sine_shaker_responses02_curr','pdf')
close(2)

%% plot RA absolute and tuning thresholds

figure(3)
loglog(freqs(1,:),thresRA,'k','LineWidth',2)
hold on
loglog(freqs(1,:),tuningRA,'r','LineWidth',2)
title('RA thresholds')
legend('Absolute threshold','Tuning threshold')

screenshot(3,'figs/sine_shaker_responses03_curr','pdf')
close(3)

%% plot RA tuning thresholds

figure(4)
loglog(freqs(1:9,:)',tuning(a.iRA,:)','k','LineWidth',2)
title('Tuning curves for RAs')
xlim([1 500])

screenshot(4,'figs/sine_shaker_responses04_curr','pdf')
close(4)

%% plot PC absolute and tuning thresholds

figure(5)
loglog(freqs(1,:),thresPC,'k','LineWidth',2)
hold on
loglog(freqs(1,:),tuningPC,'r','LineWidth',2)
title('PC thresholds')
legend('Absolute threshold','Tuning threshold')

screenshot(5,'figs/sine_shaker_responses05_curr','pdf')
close(5)

%% plot PC tuning thresholds

figure(6)
loglog(freqs(1:4,:)',tuning(a.iPC,:)','k','LineWidth',2)
title('Tuning curves for PCs')

screenshot(6,'figs/sine_shaker_responses06_curr','pdf')
close(6)
