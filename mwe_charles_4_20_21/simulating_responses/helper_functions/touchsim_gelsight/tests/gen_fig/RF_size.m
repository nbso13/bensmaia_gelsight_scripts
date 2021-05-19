function RF_size()

numSA = 4;
numRA = 9;
numPC = 4;

%% SA1

depths = linspace(0.005,0.3,25);

for n=1:numSA
    min_loc = 0;
    max_loc = 30;
    
    % find threshold
    a = Afferent('SA1','idx',n);
    for d=1:length(depths)
        s = stim_probe(depths(d),40);
        %s = stim_ramp(depths(d),0.15,[],[],[],'sine',0.25);
        r = a.response(s);
        rr_thres(d) = r.rate;
    end
    thres_ind = find(rr_thres>1,1);
    s = stim_probe(depths(thres_ind)*5,40);
    %s = stim_ramp(depths(thres_ind),0.15,[],[],[],'sine',0.25);
    
    for i=1:10
        a = affpop_linear((max_loc-min_loc)/2,max_loc-min_loc,'SA1',n);
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
    RFrad_SA(n) = (max_loc+min_loc)/2;
end

numSA = length(RFrad_SA);

%% RA

for n=1:numRA
    min_loc = 0;
    max_loc = 50;
    
    % find threshold
    a = Afferent('RA','idx',n);
    for d=1:length(depths)
        s = stim_probe(depths(d),40);
        r = a.response(s);
        rr_thres(d) = r.rate;
    end
    thres_ind = find(rr_thres>1,1);
    s = stim_probe(depths(thres_ind)*5,40);
    
    for i=1:10
        a = affpop_linear((max_loc-min_loc)/2,max_loc-min_loc,'RA',n);
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

numRA = length(RFrad_RA);

%% PC

s = stim_probe(0.05,250);

for n=1:numPC
    min_loc = 0;
    max_loc = 150;
    
    for i=1:10
        a = affpop_linear((max_loc-min_loc)/2,max_loc-min_loc,'PC',n);
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

numPC = length(RFrad_PC);

%% generate RF size figure

a = affpop_hand();

close all
figure(1)
set(gcf, 'Position', [300    1   300   1080]);

ax = subplot(311);
[origin,theta,pxl_per_mm,ext_space] = plot_hand(ax,'names',0,'axes',0,'centers',0);
rot=[cos(-theta) -sin(-theta);sin(-theta) cos(-theta)];
indSA = randperm(length(a.iSA1));
locSA = a.location(indSA(1:numSA),:)*rot*pxl_per_mm + repmat(origin,numSA,1);
viscircles(locSA,RFrad_SA+eps);
title('SA1')
axis equal

ax = subplot(312);
plot_hand(ax,'names',0,'axes',0,'centers',0);
indRA = randperm(length(a.iRA));
locRA = a.location(indRA(1:numRA),:)*rot*pxl_per_mm + repmat(origin,numRA,1);
idx_ok = ~isnan(RFrad_RA);
viscircles(locRA(idx_ok,:),RFrad_RA(idx_ok)+eps);
plot(locRA(~idx_ok,1),locRA(~idx_ok,1),'rx','MarkerSize',15)
title('RA')
axis equal

ax = subplot(313);
plot_hand(ax,'names',0,'axes',0,'centers',0);
indPC = randperm(length(a.iPC));
locPC = a.location(indPC(1:numPC),:)*rot*pxl_per_mm + repmat(origin,numPC,1);
viscircles(locPC,RFrad_PC+eps);
title('PC')
axis equal


screenshot(1,'figs/RF_size01_curr','pdf')
close(1)
clear r

%% growth of RF with indentation depth

indent_depths = [50 100 200 350 500]/1000;

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



figure(2)
set(gcf,'pos',[0 0 480 948]);
subplot(211)
hold on
plot(indent_depths,RFsiz_SA,'Color',[0 .5 0],'LineWidth',2)
plot(indent_depths,RFsiz_RA,'Color',[0 0 .5],'LineWidth',2)
subplot(212)
hold on
plot(indent_depths,nanmean(RFsiz_SA,2),'-o','Color',[0 .5 0],'LineWidth',2)
plot(indent_depths,nanmean(RFsiz_RA,2),'-o','Color',[0 0 .5],'LineWidth',2)
legend('SA','RA','Location','NorthWest')

screenshot(2,'figs/RF_size02_curr','pdf')
close(2)
