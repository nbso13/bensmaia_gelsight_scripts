% receptive_fields.m

cols=affcol;

%% determine RF size for the different afferent classes

numSA = 4;
numRA = 9;
numPC = 4;

% SA1
depths = linspace(0.005,0.3,25);
for n=1:numSA
    min_loc = 0;
    max_loc = 30;
    
    % find threshold
    a = Afferent('SA1','idx',n,'noisy',false);
    for d=1:length(depths)
        s = stim_probe(depths(d),25);
        %s = stim_ramp(depths(d),0.15,[],[],0.005,'sine',0.5);
        r = a.response(s);
        rr_thres(d) = r.rate;
    end
    thres_ind = find(rr_thres>1,1);
    s = stim_probe(depths(thres_ind)*5,25);
    %s = stim_ramp(depths(thres_ind)*5,0.15,[],[],0.005,'sine',0.5);
    
    for i=1:10
        a = affpop_linear((max_loc-min_loc)/2,max_loc-min_loc,'SA1',n,'noisy',false);
        a.location = a.location + repmat([min_loc 0],3,1);
        for rep=1:5
            r = a.response(s);
            rr(:,rep) = r.rate;
        end
        idx = find(mean(rr,2)<=8,1);
        if isempty(idx) || idx==1
            break;
        end
        max_loc = a.location(idx,1);
        min_loc = a.location(idx-1,1);
    end
    RFrad_SA(n) = (max_loc+min_loc)/2+0.25;
end

% RA
for n=1:numRA
    min_loc = 0;
    max_loc = 50;
    
    % find threshold
    a = Afferent('RA','idx',n,'noisy',false);
    for d=1:length(depths)
        s = stim_probe(depths(d),25);
        %s = stim_ramp(depths(d),0.15,[],[],0.01,'sine',0.5);
        r = a.response(s);
        rr_thres(d) = r.rate;
    end
    thres_ind = find(rr_thres>1,1);
    s = stim_probe(depths(thres_ind)*5,25);
    %s = stim_ramp(depths(thres_ind)*5,0.15,[],[],0.01,'sine',0.5);
    
    for i=1:10
        a = affpop_linear((max_loc-min_loc)/2,max_loc-min_loc,'RA',n,'noisy',false);
        a.location = a.location + repmat([min_loc 0],3,1);
        for rep=1:5
            r = a.response(s);
            rr(:,rep) = r.rate;
        end
        idx = find(mean(rr,2)<=8,1);
        if isempty(idx) || idx==1
            break;
        end
        max_loc = a.location(idx,1);
        min_loc = a.location(idx-1,1);
    end
    RFrad_RA(n) = (max_loc+min_loc)/2;
end

% PC
depths = linspace(0.0001,0.025,25);
for n=1:numPC
    min_loc = 0;
    max_loc = 150;
    
    % find threshold
    a = Afferent('PC','idx',n,'noisy',false);
    for d=1:length(depths)
        s = stim_probe(depths(d),250);
        %s = stim_ramp(depths(d),0.15,[],[],0.005,'sine',0.5);
        r = a.response(s);
        rr_thres(d) = r.rate;
    end
    thres_ind = find(rr_thres>1,1);
    s = stim_probe(depths(thres_ind)*5,250);
    %s = stim_ramp(depths(thres_ind)*5,0.15,[],[],0.005,'sine',0.5);
    
    for i=1:10
        a = affpop_linear((max_loc-min_loc)/2,max_loc-min_loc,'PC',n,'noisy',false);
        a.location = a.location + repmat([min_loc 0],3,1);
        for rep=1:5
            r = a.response(s);
            rr(:,rep) = r.rate;
        end
        idx = find(r.rate<=8,1);
        if isempty(idx) || idx==1
            break;
        end
        max_loc = a.location(idx,1);
        min_loc = a.location(idx-1,1);
    end
    RFrad_PC(n) = (max_loc+min_loc)/2;
end

%% generate RF size figure

%source: Vallbo & Johansson 1980
XSA = csvread('receptive_fields1a.csv',1);
XRA = csvread('receptive_fields1b.csv',1);
XPC = csvread('receptive_fields1c.csv',1);
pcSA = interp1(XSA(:,2),XSA(:,1),[.25 .5 .75]);
pcRA = interp1(XRA(:,2),XRA(:,1),[.25 .5 .75]);
pcPC = interp1(XPC(:,2),XPC(:,1),[.25 .5 .75]);
pcSA(2) = 11; pcRA(2) = 12.6; pcPC(2) = 101.3;

figure(1)
clf
set(gcf,'pos',[0 0 240 474]);

subplot(211)

line([0.6 1.4],[pcSA(1) pcSA(1)],'LineWidth',1.5,'Color',affcol(1))
line([0.6 1.4],[pcSA(3) pcSA(3)],'LineWidth',1.5,'Color',affcol(1))
line([1 1],[pcSA(1) pcSA(3)],'LineWidth',1.5,'Color',affcol(1))
line([0.5 1.5],[pcSA(2) pcSA(2)],'LineWidth',1.5,'Color','k')

line([1.6 2.4],[pcRA(1) pcRA(1)],'LineWidth',1.5,'Color',affcol(2))
line([1.6 2.4],[pcRA(3) pcRA(3)],'LineWidth',1.5,'Color',affcol(2))
line([2 2],[pcRA(1) pcRA(3)],'LineWidth',1.5,'Color',affcol(2))
line([1.5 2.5],[pcRA(2) pcRA(2)],'LineWidth',1.5,'Color','k')

line([2.6 3.4],[pcPC(1) pcPC(1)],'LineWidth',1.5,'Color',affcol(3))
line([2.6 3.4],[pcPC(3) pcPC(3)],'LineWidth',1.5,'Color',affcol(3))
line([3 3],[pcPC(1) pcPC(3)],'LineWidth',1.5,'Color',affcol(3))
line([2.5 3.5],[pcPC(2) pcPC(2)],'LineWidth',1.5,'Color','k')
set(gca,'yscale','log','ylim',[1 550],'xtick',[1 2 3])
ylabel('RF area [mm^2]')
box off

subplot(212)
semilogy(1+randn(1,numSA)/5,pi*RFrad_SA.^2,'o','Color',cols(1,:),'LineWidth',2)
hold on
line([0.5 1.5],[mean(pi*RFrad_SA.^2) mean(pi*RFrad_SA.^2)],'LineWidth',1.5,'Color','k')
semilogy(2+randn(1,numRA)/5,pi*RFrad_RA.^2,'o','Color',cols(2,:),'LineWidth',2)
line([1.5 2.5],[mean(pi*RFrad_RA.^2) mean(pi*RFrad_RA.^2)],'LineWidth',1.5,'Color','k')
semilogy(3+randn(1,numPC)/5,pi*RFrad_PC.^2,'o','Color',cols(3,:),'LineWidth',2)
line([2.5 3.5],[mean(pi*RFrad_PC.^2) mean(pi*RFrad_PC.^2)],'LineWidth',1.5,'Color','k')
set(gca,'yscale','log','ylim',[1 550],'xtick',[1 2 3])
ylabel('RF area [mm^2]')
box off

screenshot(1,'receptive_fields1','pdf')
close(1)

%% growth of RF with indentation depth for SAs and RAs

indent_depths = [50 100 200 350 500]/1000;
r = [];

for d=1:length(indent_depths)
    s = stim_ramp(indent_depths(d),0.2,[],[],0.01,'sine',0.25,1.7);
    
    for n=1:numSA
        a = affpop_linear(.5,5,'SA1',n,'noisy',true);
        for rep=1:5
            r(:,rep) = a.response(s).rate;
        end
        r = mean(r,2);
        if r(1)==0
            RFsiz_SA(d,n) = NaN;
            continue;
        end        
        num_probes = sum(r>=0.1*r(1));
        RFsiz_SA(d,n) = (num_probes*0.5)^2*pi;
    end
    
    for n=1:numRA
        a = affpop_linear(.5,5,'RA',n,'noisy',true);
        r = a.response(s).rate;
        if r(1)==0
            RFsiz_RA(d,n) = NaN;
            continue;
        end        
        num_probes = sum(r>=0.1*r(1));
        RFsiz_RA(d,n) = (num_probes*0.5)^2*pi;
    end
end

RFsiz_SA(:,isnan(sum(RFsiz_SA))) = [];
RFsiz_RA(:,isnan(sum(RFsiz_RA))) = [];

% read in published numbers
X = csvread('receptive_fields2.csv',1);

figure(2)
set(gcf,'pos',[0 0 240 474]);
subplot(211)
hold on
plot(X(1:2:end,1)/1000,X(1:2:end,2),'-o','Color',cols(1,:),'LineWidth',2)
plot(X(1:2:end,1)/1000,X(1:2:end,3),'-o','Color',cols(2,:),'LineWidth',2)
ylim([0 25])
legend('SA1','RA','Location','NorthWest')
ylabel('RF area [mm^2]')
subplot(212)
hold on
plot(indent_depths,nanmean(RFsiz_SA,2),'-o','Color',cols(1,:),'LineWidth',2)
plot(indent_depths,nanmean(RFsiz_RA,2),'-o','Color',cols(2,:),'LineWidth',2)
xlabel('Indentation depth [um]')
ylabel('RF area [mm^2]')
ylim([0 25])

screenshot(2,'receptive_fields2','pdf')
close(2)

%% amplitude needed to generate response away from RF center for RAs

X = csvread('receptive_fields3.csv',1);

%locs = logspace(-1,1,10);
locs = X(:,1)';
nlocs=length(locs);

ap = AfferentPopulation;
for l=1:nlocs
    a_sub = affpop_single_models([locs(l) 0]);
    a_sub.afferents(~a_sub.iRA)=[]; % only RA's
    ap.afferents = [ap.afferents a_sub.afferents];
end
naff=length(ap.afferents);

% create stimuli and compute resp
frq=40;
rad=1;
pre_indent=.55;
amp = linspace(0,1,60);
namp=length(amp);
resp_rate=zeros(namp,naff);
for ii=1:namp
    stim=stim_sine(frq,amp(ii),[],[],[],[],[],rad,pre_indent);
    rc(ii)=ap.response(stim);
    resp_rate(ii,:)=rc(ii).rate;
    fprintf('.')
end
fprintf('\n')


idx_phaseloc=zeros(naff,1);
for ii=1:naff
    % find minimal (first) ampl showing phase locked resp (.9 spikes/cycle)
    truc=find(resp_rate(:,ii)>frq*.9,1);
    if(isempty(truc)),       truc=namp;   end
    idx_phaseloc(ii)=truc;
end

% reshape phaseloc resp to matrix nmodels x nlocs
amp_phaseloc=reshape(amp(idx_phaseloc),naff/nlocs,nlocs);

% Fig 7 Johnson 1974
figure(3)
set(gcf,'pos',[0 0 240 474]);
subplot(211)
loglog(X(:,1),X(:,2),'ok','LineWidth',2,'Color',cols(2,:))
set(gca,'xscale','log','yscale','log','xlim',[.1 10],'ylim',[10 1000])
ylabel('Minimum amplitude [um]')
box off

subplot(212)
loglog(locs,nanmean(amp_phaseloc)*1000,'ok','LineWidth',2,'Color',cols(2,:))
set(gca,'xscale','log','yscale','log','xlim',[.1 10],'ylim',[10 1000])
xlabel('Distance from stimulus [mm]')
ylabel('Minimum amplitude [um]')
box off

screenshot(3,'receptive_fields3','pdf')
close(3)

%% RF sizes for high-amplitude stimuli

numSA = 4;
numRA = 9;
numPC = 4;

% SA1
for n=1:numSA
    min_loc = 0;
    max_loc = 30;
    
    % find threshold
    a = Afferent('SA1','idx',n,'noisy',false);
    s = stim_probe(0.75,25);
    %s = stim_ramp(depths(thres_ind)*5,0.15,[],[],0.005,'sine',0.5);
    
    for i=1:10
        a = affpop_linear((max_loc-min_loc)/2,max_loc-min_loc,'SA1',n,'noisy',false);
        a.location = a.location + repmat([min_loc 0],3,1);
        for rep=1:5
            r = a.response(s);
            rr(:,rep) = r.rate;
        end
        idx = find(mean(rr,2)<=8,1);
        if isempty(idx) || idx==1
            break;
        end
        max_loc = a.location(idx,1);
        min_loc = a.location(idx-1,1);
    end
    RFradHI_SA(n) = (max_loc+min_loc)/2+0.25;
end

% RA
for n=1:numRA
    min_loc = 0;
    max_loc = 50;
    
    % find threshold
    a = Afferent('RA','idx',n,'noisy',false);
    s = stim_probe(0.75,25);
    %s = stim_ramp(depths(thres_ind)*5,0.15,[],[],0.01,'sine',0.5);
    
    for i=1:10
        a = affpop_linear((max_loc-min_loc)/2,max_loc-min_loc,'RA',n,'noisy',false);
        a.location = a.location + repmat([min_loc 0],3,1);
        for rep=1:5
            r = a.response(s);
            rr(:,rep) = r.rate;
        end
        idx = find(mean(rr,2)<=8,1);
        if isempty(idx) || idx==1
            break;
        end
        max_loc = a.location(idx,1);
        min_loc = a.location(idx-1,1);
    end
    RFradHI_RA(n) = (max_loc+min_loc)/2;
end

% PC
for n=1:numPC
    min_loc = 0;
    max_loc = 150;
    
    % find threshold
    a = Afferent('PC','idx',n,'noisy',false);
    s = stim_probe(0.03,250);
    %s = stim_ramp(depths(thres_ind)*5,0.15,[],[],0.005,'sine',0.5);
    
    for i=1:10
        a = affpop_linear((max_loc-min_loc)/2,max_loc-min_loc,'PC',n,'noisy',false);
        a.location = a.location + repmat([min_loc 0],3,1);
        for rep=1:5
            r = a.response(s);
            rr(:,rep) = r.rate;
        end
        idx = find(r.rate<=8,1);
        if isempty(idx) || idx==1
            break;
        end
        max_loc = a.location(idx,1);
        min_loc = a.location(idx-1,1);
    end
    RFradHI_PC(n) = (max_loc+min_loc)/2;
end

% generate RF size figure
figure(4)
set(gcf,'pos',[0 0 240 474]);
subplot(211)
loglog(pi*RFrad_SA.^2,pi*RFradHI_SA.^2,'o','Color',cols(1,:))
hold on
loglog(pi*RFrad_RA.^2,pi*RFradHI_RA.^2,'o','Color',cols(2,:))
loglog(pi*RFrad_PC.^2,pi*RFradHI_PC.^2,'o','Color',cols(3,:))
unityslope;

subplot(212)
semilogy(1+randn(1,numSA)/5,pi*RFradHI_SA.^2,'o','Color',cols(1,:))
hold on
line([0.5 1.5],[mean(pi*RFradHI_SA.^2) mean(pi*RFradHI_SA.^2)],'LineWidth',1.5,'Color','k')
semilogy(2+randn(1,numRA)/5,pi*RFradHI_RA.^2,'o','Color',cols(2,:))
line([1.5 2.5],[mean(pi*RFradHI_RA.^2) mean(pi*RFradHI_RA.^2)],'LineWidth',1.5,'Color','k')
semilogy(3+randn(1,numPC)/5,pi*RFradHI_PC.^2,'o','Color',cols(3,:))
line([2.5 3.5],[mean(pi*RFradHI_PC.^2) mean(pi*RFradHI_PC.^2)],'LineWidth',1.5,'Color','k')
set(gca,'yscale','log','ylim',[1 10e5])
ylabel('RF area [mm^2]')
box off

screenshot(4,'receptive_fields4','pdf')
close(4)

%% show PCs responding to stimuli far away

a = affpop_hand();

% 250 Hz stim
s1 = stim_sine(250,0.015,[],0.1);
r1 = a.response(s1,1);
ind_resp1PC = r1.rate(a.iPC)>0;
ind_resp1RA = r1.rate(a.iRA)>0;

% 30 Hz stim
s2 = stim_sine(30,0.25,[],0.1);
r2 = a.response(s2,1);
ind_resp2PC = r2.rate(a.iPC)>0;
ind_resp2RA = r2.rate(a.iRA)>0;
%%

figure(5)
set(gcf,'pos',[0 0 240 474]);
subplot(211)
[origin,theta,pxl_per_mm] = plot_hand(gca,'names',false,'axes',false,'centers',false);

all_locsPC = a.location(a.iPC,:);
rot = [cos(-theta) -sin(-theta);sin(-theta) cos(-theta)];
all_locsPC = all_locsPC*rot;
all_locsPC = all_locsPC*pxl_per_mm;
all_locsPC = bsxfun(@plus,all_locsPC,origin);

all_locsRA = a.location(a.iRA,:);
all_locsRA = all_locsRA*rot;
all_locsRA = all_locsRA*pxl_per_mm;
all_locsRA = bsxfun(@plus,all_locsRA,origin);

plot(all_locsPC(:,1),all_locsPC(:,2),'.','Color',[.75 .75 .75])
plot(all_locsPC(ind_resp1PC,1),all_locsPC(ind_resp1PC,2),'.','Color',cols(3,:))


subplot(212)
plot_hand(gca,'names',false,'axes',false,'centers',false);

plot(all_locsRA(:,1),all_locsRA(:,2),'.','Color',[.75 .75 .75])
plot(all_locsRA(ind_resp2RA,1),all_locsRA(ind_resp2RA,2),'.','Color',cols(2,:))

screenshot(5,'receptive_fields5','pdf')
close(5)


figure(6)
set(gcf,'pos',[0 0 240 474]);

subplot(211)
plot(s1.trace)
axis off

subplot(212)
plot(s2.trace)
axis off

screenshot(6,'receptive_fields6','pdf')
close(6)
