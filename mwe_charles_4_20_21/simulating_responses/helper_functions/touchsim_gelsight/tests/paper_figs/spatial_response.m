% spatial_response.m
cols=affcol;

%% show firing rate versus number of probes

% set up afferent population
aSA = affpop_single_models([0 0],'SA1');
aRA = affpop_single_models([0 0],'RA');

% set up stimulus and generate responses
cpin = [0 0];
pins = [ 1 -0.5;
         1  0.5;
         0 -1;
         0  1;
        -1 -0.5;
        -1  0.5];

s{1} = stim_ramp(0.1,0.2,cpin,[],0.0125,[],0.25,1.5);
rSA(:,1) = aSA.response(s{1}).rate;
rRA(:,1) = aRA.response(s{1}).rate;
for n=1:6
    perms = nchoosek(1:6,n);
    rSA_tmp = [];
    rRA_tmp = [];
    for i=1:size(perms,1)
        s{n+1}(i) = stim_indent_shape([cpin; pins(perms(i,:),:)],stim_ramp(0.1,0.2,cpin,[],0.0125,[],0.25,1.5));
        rSA_tmp(:,i) = aSA.response(s{n+1}(i)).rate;
        rRA_tmp(:,i) = aRA.response(s{n+1}(i)).rate;
    end
    rSA(:,n+1) = mean(rSA_tmp,2);
    rRA(:,n+1) = mean(rRA_tmp,2);
end

% SA plots

X = csvread('spatial_response1.csv',1);

figure(1)
set(1,'pos',[0 0 240 474])

subplot(211)
semilogy(1:7,X(:,2),'-o','LineWidth',1.5,'Color',cols(1,:))
ylim([1 500])
xlim([0 7])
box off
ylabel('Firing rate [Hz]')

subplot(212)
semilogy(1:7,mean(rSA),'-o','LineWidth',1.5,'Color',cols(1,:))
ylim([1 500])
xlim([0 7])
box off
xlabel('Number of probes')
ylabel('Firing rate [Hz]')

screenshot(1,'spatial_response1','pdf')
close(1)


% RA plots

X = csvread('spatial_response2.csv',1);

figure(2)
set(2,'pos',[0 0 240 474])

subplot(211)
semilogy(1:7,X(:,2),'-o','LineWidth',1.5,'Color',cols(2,:))
ylim([1 500])
xlim([0 7])
box off
ylabel('Firing rate [Hz]')

subplot(212)
semilogy(1:7,mean(rRA),'-ok','LineWidth',1.5,'Color',cols(2,:))
ylim([1 500])
xlim([0 7])
box off
xlabel('Number of probes')
ylabel('Firing rate [Hz]')

screenshot(2,'spatial_response2','pdf')
close(2)


%% show edge enhancement for different bars

img = [nan(1,15) zeros(1,30) nan(1,5) zeros(1,5) nan(1,5) zeros(1,5) nan(1,8) zeros(1,5) nan(1,10) zeros(1,5) nan(1,15) zeros(1,5) nan(1,20) zeros(1,5) nan(1,30) zeros(1,5) nan(1,50) zeros(1,30) nan(1,15)]';
img= repmat(img,1,30);
shape = img2shape(img,[0 NaN],10);
shape(:,1)=shape(:,1)+1.5;

s = stim_indent_shape(shape,stim_ramp(2));

% compute response
ap=affpop_linear(0.2,28.4);           % create full affpop
ap.afferents(ap.iPC)=[];            % remove PCs
modelidx=[ap.afferents(:).idx]';    % get model indices
rc=ap.response(s);                  % compute response

% parse rates
allrates=rc.rate;
% loop over SA1s
for n=1:length(unique(modelidx(ap.iSA1)))
    rates_SA(:,n) = allrates(modelidx==n&ap.iSA1');
end
% loop over RAs
for n=1:length(unique(modelidx(ap.iRA)))
    rates_RA(:,n) = allrates(modelidx==n&ap.iRA');
end

% SA plot
X1 = csvread('spatial_response1a.csv',1);
X2 = csvread('spatial_response1b.csv',1);
X3 = csvread('spatial_response1c.csv',1);
X4 = csvread('spatial_response1d.csv',1);

stimxz=(0:285)'/10;
[~,i]=intersect(stimxz(:,1),shape(:,1));
stimxz(i,2)=1;

figure(3)
set(gcf,'pos',[0 0 480 474]);
subplot(2,1,1), hold on
plot(stimxz(:,1),stimxz(:,2)*15,'-','color',[0 0 0 .7],'LineWidth',1.5)
plot(X1(:,1),X1(:,2),'LineWidth',1.5,'Color',cols(1,:))
plot(X2(:,1),X2(:,2),'LineWidth',1.5,'Color',cols(1,:))
plot(X3(:,1),X3(:,2),'LineWidth',1.5,'Color',cols(1,:))
plot(X4(:,1),X4(:,2),'LineWidth',1.5,'Color',cols(1,:))
xlim([0 27])
ylim([0 100])
ylabel('Firing rate [Hz]')

subplot(2,1,2), hold on
plot(stimxz(:,1),stimxz(:,2)*15,'-','color',[0 0 0 .7],'LineWidth',1.5)
for n=1:length(unique(modelidx(ap.iSA1)))
    idx=modelidx==n&ap.iSA1';
    plot(ap.location(idx,1),rates_SA(:,n),'LineWidth',1.5,'Color',cols(1,:))
end
xlim([0 27])
ylim([0 100])
ylabel('Firing rate [Hz]')
xlabel('Position [mm]')

screenshot(3,'spatial_response3','pdf')
close(3)


% RA plot
X1 = csvread('spatial_response2a.csv',1);
X2 = csvread('spatial_response2b.csv',1);
X3 = csvread('spatial_response2c.csv',1);
X4 = csvread('spatial_response2d.csv',1);
X5 = csvread('spatial_response2e.csv',1);

figure(4)
set(gcf,'pos',[0 0 480 474]);
subplot(2,1,1), hold on
plot(stimxz(:,1),stimxz(:,2)*3,'-','color',[0 0 0 .7],'LineWidth',1.5)
plot(X1(:,1),X1(:,2),'LineWidth',1.5,'Color',cols(2,:))
plot(X2(:,1),X2(:,2),'LineWidth',1.5,'Color',cols(2,:))
plot(X3(:,1),X3(:,2),'LineWidth',1.5,'Color',cols(2,:))
plot(X4(:,1),X4(:,2),'LineWidth',1.5,'Color',cols(2,:))
plot(X5(:,1),X5(:,2),'LineWidth',1.5,'Color',cols(2,:))
xlim([0 27])
ylim([0 20])
ylabel('Firing rate [Hz]')

subplot(2,1,2), hold on
plot(stimxz(:,1),stimxz(:,2)*3,'-','color',[0 0 0 .7],'LineWidth',1.5)
for n=1:length(unique(modelidx(ap.iSA1)))
    idx=modelidx==n&ap.iRA';
    plot(ap.location(idx,1),rates_RA(:,n),'LineWidth',1.5,'Color',cols(2,:))
end
xlim([0 27])
ylim([0 20])
ylabel('Firing rate [Hz]')
xlabel('Position [mm]')

screenshot(4,'spatial_response4','pdf')
close(4)
