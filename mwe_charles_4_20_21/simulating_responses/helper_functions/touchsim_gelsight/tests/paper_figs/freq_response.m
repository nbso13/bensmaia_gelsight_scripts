% freq_response.m

cols=affcol;


%% set up stimuli
freqs = [1 5 10 25 50 100 200 300 400 600 800 1000];
freqs = repmat(freqs,14,1);

amps(:,1) = [6.9477    9.6541   13.4148   18.6405   25.9017   35.9915   50.0118   69.4934   96.5641  134.1799  186.4486  259.0781  360.0000  500.0000]';
amps(:,2) = [6.9477    9.6541   13.4148   18.6405   25.9017   35.9915   50.0118   69.4934   96.5641  134.1799  186.4486  259.0781  360.0000  500.0000]';
amps(:,3) = [3.3880    4.5914    6.2222    8.4322   11.4272   15.4861   20.9866   28.4408   38.5427   52.2327   70.7853   95.9275  130.0000  259.0781]';
amps(:,4) = [3.3880    4.5914    6.2222    8.4322   11.4272   15.4861   20.9866   28.4408   38.5427   52.2327   70.7853   95.9275  130.0000  259.0781]';
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

thres_num_spikes = 1;

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

%% compare parameters of mean rate-intensity functions

figure(1)
set(gcf,'pos',[0 0 240 474]);

X1 = csvread('freq_response1a.csv',1);
X2 = csvread('freq_response1b.csv',1);
X3 = csvread('freq_response1c.csv',1);

subplot(211)
loglog(X1(:,1),X1(:,2),'LineWidth',2,'Color',cols(1,:))
hold on
loglog(X2(:,1),X2(:,2),'LineWidth',2,'Color',cols(2,:))
loglog(X3(:,1),X3(:,2),'LineWidth',2,'Color',cols(3,:))
ylim([1e0 1e3])
xlim([1 1000])
set(gca,'xtick',[1 10 100 1000])
legend('SA1','RA','PC','Location','NorthWest')
ylabel('\alpha')
box off

subplot(212)
loglog(freqs(1,:),pSA(:,1),'Color',cols(1,:),'LineWidth',2);
hold on
loglog(freqs(1,:),pRA(:,1),'Color',cols(2,:),'LineWidth',2);
loglog(freqs(1,:),pPC(:,1),'Color',cols(3,:),'LineWidth',2);
ylim([1e0 1e3])
xlim([1 1000])
set(gca,'xtick',[1 10 100 1000])
ylabel('\alpha')
box off

screenshot(1,'freq_response1','pdf')
close(1)


figure(2)
set(gcf,'pos',[0 0 240 474]);

X1 = csvread('freq_response2a.csv',1);
X2 = csvread('freq_response2b.csv',1);
X3 = csvread('freq_response2c.csv',1);

subplot(211)
loglog(X1(:,1),X1(:,2),'LineWidth',2,'Color',cols(1,:))
hold on
loglog(X2(:,1),X2(:,2),'LineWidth',2,'Color',cols(2,:))
loglog(X3(:,1),X3(:,2),'LineWidth',2,'Color',cols(3,:))
ylim([1e-1 1e2])
xlim([1 1000])
set(gca,'xtick',[1 10 100 1000])
legend('SA1','RA','PC','Location','NorthWest')
ylabel('\alpha')
box off

subplot(212)
loglog(freqs(1,:),10.^pSA(:,2),'Color',cols(1,:),'LineWidth',2);
hold on
loglog(freqs(1,:),10.^pRA(:,2),'Color',cols(2,:),'LineWidth',2);
loglog(freqs(1,:),10.^pPC(:,2),'Color',cols(3,:),'LineWidth',2);
ylim([1e-1 1e2])
xlim([1 1000])
set(gca,'xtick',[1 10 100 1000])
xlabel('Frequency (Hz)')
ylabel('10^\beta (\mum)')
box off

screenshot(2,'freq_response2','pdf')
close(2)

%% plot RA absolute and tuning thresholds

figure(3)

set(gcf,'pos',[0 0 240 474]);
subplot(211)

subplot(212)
loglog(freqs(1,:),thresRA,'k','LineWidth',2)
hold on
loglog(freqs(1,:),tuningRA,'r','LineWidth',2)
legend('Absolute thres.','Tuning thres.')
xlabel('Frequency [Hz]')
ylabel('Amplitude [um]')
box off

screenshot(3,'freq_response3','pdf')
close(3)

%% plot RA tuning thresholds

X = csvread('freq_response4.csv',1);

for n=2:12
    idx = find(X(:,n)==X(1,n));
    X(idx(1:end-1),n) = NaN;
    idx = find(X(:,n)==X(end,n));
    X(idx(2:end),n) = NaN;
end

figure(4)
set(gcf,'pos',[0 0 240 474]);
subplot(211)
loglog(X(:,1),X(:,2:end),'k','LineWidth',2)
xlim([1 500])
ylabel('Amplitude [um]')
box off

subplot(212)
loglog(freqs(1:9,:)',tuning(a.iRA,:)','k','LineWidth',2)
xlim([1 500])
xlabel('Frequency [Hz]')
ylabel('Amplitude [um]')
box off

screenshot(4,'freq_response4','pdf')
close(4)

%% plot PC absolute and tuning thresholds

figure(5)
set(gcf,'pos',[0 0 240 474]);
subplot(211)

subplot(212)
loglog(freqs(1,:),thresPC,'k','LineWidth',2)
hold on
loglog(freqs(1,:),tuningPC,'r','LineWidth',2)
legend('Absolute threshold','Tuning threshold')
xlabel('Frequency [Hz]')
ylabel('Amplitude [um]')
box off

screenshot(5,'freq_response5','pdf')
close(5)

%% plot PC tuning thresholds

X = csvread('freq_response6.csv',1);

for n=2:12
    idx = find(X(:,n)==X(1,n));
    X(idx(1:end-1),n) = NaN;
    idx = find(X(:,n)==X(end,n));
    X(idx(2:end),n) = NaN;
end

figure(6)
set(gcf,'pos',[0 0 240 474]);
subplot(211)
loglog(X(:,1),X(:,2:end),'k','LineWidth',2)
xlim([10 500])
ylim([.5 1000])
ylabel('Amplitude [um]')
box off

subplot(212)
loglog(freqs(1:4,:)',tuning(a.iPC,:)','k','LineWidth',2)
xlim([10 500])
ylim([.5 1000])
xlabel('Frequency [Hz]')
ylabel('Amplitude [um]')
box off

screenshot(6,'freq_response6','pdf')
close(6)

%%

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

%% plot johansson freq response figs
SAcol=[199,233,192;161,217,155;116,196,118;65,171,93;35,139,69;0,109,44;0,68,27]/255;
RAcol=[198,219,239;158,202,225;107,174,214;66,146,198;33,113,181;8,81,156;8,48,107]/255;
PCcol=[253,208,162;253,174,107;253,141,60;241,105,19;217,72,1;166,54,3;127,39,4]/255;

spkpercycleSA=median(spkpercycle(:,:,ap.iSA1),3);
spkpercycleRA=median(spkpercycle(:,:,ap.iRA),3);
spkpercyclePC=median(spkpercycle(:,:,ap.iPC),3);

% SA1
X = csvread('freq_response7.csv',1);
for n=2:8
    idx = find(X(:,n)==X(1,n));
    X(idx(1:end-1),n) = NaN;
    idx = find(X(:,n)==X(end,n));
    X(idx(2:end),n) = NaN;
end


figure(7); set(0,'defaultAxesColorOrder',flipud(SAcol))
set(gcf,'pos',[0 0 240 474]);
subplot(211)
semilogx(X(:,1),X(:,2:end),'LineWidth',1.5);
set(gca,'xlim',[.4 500],'ylim',[0 7],'xtick',[1 10 100])
ylabel('Spikes per cycle')
box off

subplot(212)
semilogx(frq,spkpercycleSA,'LineWidth',1.5);
set(gca,'xlim',[.4 500],'ylim',[0 7],'xtick',[1 10 100])
xlabel('Frequency [Hz]')
ylabel('Spikes per cycle')
box off

screenshot(7,'freq_response7','pdf')
%close(7)


% RA
X = csvread('freq_response8.csv',1);
for n=2:8
    idx = find(X(:,n)==X(1,n));
    X(idx(1:end-1),n) = NaN;
    idx = find(X(:,n)==X(end,n));
    X(idx(2:end),n) = NaN;
end

figure(8); set(0,'defaultAxesColorOrder',flipud(RAcol))
set(gcf,'pos',[0 0 240 474]);
subplot(211)
semilogx(X(:,1),X(:,2:end),'LineWidth',1.5);
set(gca,'xlim',[.4 500],'ylim',[0 7],'xtick',[1 10 100])
ylabel('Spikes per cycle')
box off

subplot(212)
semilogx(frq,spkpercycleRA,'LineWidth',1.5);
set(gca,'xlim',[.4 500],'ylim',[0 7],'xtick',[1 10 100])
xlabel('Frequency [Hz]')
ylabel('Spikes per cycle')
box off

screenshot(8,'freq_response8','pdf')
%close(8)


% PC
X = csvread('freq_response9.csv',1);
for n=2:8
    idx = find(X(:,n)==X(1,n));
    X(idx(1:end-1),n) = NaN;
    idx = find(X(:,n)==X(end,n));
    X(idx(2:end),n) = NaN;
end

figure(9); set(0,'defaultAxesColorOrder',flipud(PCcol))
set(gcf,'pos',[0 0 240 474]);
subplot(211)
semilogx(X(:,1),X(:,2:end),'LineWidth',1.5);
set(gca,'xlim',[.4 500],'ylim',[0 7],'xtick',[1 10 100])
ylabel('Spikes per cycle')
box off

subplot(212)
semilogx(frq,spkpercyclePC,'LineWidth',1.5);
set(gca,'xlim',[.4 500],'ylim',[0 7],'xtick',[1 10 100])
xlabel('Frequency [Hz]')
ylabel('Spikes per cycle')
box off

screenshot(9,'freq_response9','pdf')
%close(9)

set(0,'defaultAxesColorOrder',lines)
%% freeman & johnson data
load freq_response_FJdata

for ii=1:3
    figure(9+ii)
    set(gcf,'pos',[0 0 240 474]);
    
    subplot(211)
    hold on
    plot(squeeze(FJdata{ii,1}(:,1,:)),squeeze(FJdata{ii,1}(:,2,:)),'LineWidth',1.5,'color',[.6 .6 .6 .6])
    plot(squeeze(FJdata{ii,2}(:,1,:)),squeeze(FJdata{ii,2}(:,2,:)),'LineWidth',1.5,'color',[1 .6 .6 .6])
    plot(FJdata{ii,3}(:,1,1),FJdata{ii,3}(:,2,1),'LineWidth',1.5,'color',[0 0 0])
    plot(FJdata{ii,3}(:,1,2),FJdata{ii,3}(:,2,2),'LineWidth',1.5,'color',[1 0 0])
    set(gca,'xlim',[1 300],'ylim',[.1 3000],'xtick',[1 3 10 30 100 300],...
        'xscale','log','yscale','log','ytick',[.1 1 10 100 1000])
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [um]')
    box off

    subplot(212)
    switch ii
        case 1
            ii_sub = a.iSA1;
        case 2
            ii_sub = a.iRA;
        case 3
            ii_sub = a.iPC;
    end
    hold on
    plot(freqs(1,:),thres(ii_sub,:),'LineWidth',1.5,'color',[.6 .6 .6 .6])
    plot(freqs(1,:),tuning(ii_sub,:),'LineWidth',1.5,'color',[1 .6 .6 .6])
    plot(freqs(1,:),geomean(thres(ii_sub,:)),'LineWidth',1.5,'color',[0 0 0])
    plot(freqs(1,:),geomean(tuning(ii_sub,:)),'LineWidth',1.5,'color',[1 0 0])
    
    set(gca,'xlim',[1 300],'ylim',[.1 3000],'xtick',[1 3 10 30 100 300],...
        'xscale','log','yscale','log','ytick',[.1 1 10 100 1000])
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [um]')
    box off
    
    screenshot(9+ii,['freq_response' num2str(9+ii)],'pdf')
    close(9+ii)
end

