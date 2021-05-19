function spatial_res_bars()
%clear;clc; close all
% set up bars
%img = [nan(1,5) zeros(1,5) nan(1,1) zeros(1,1) nan(1,1) zeros(1,1) nan(1,2) zeros(1,1) nan(1,3) zeros(1,1) nan(1,4) zeros(1,1) nan(1,6) zeros(1,1) nan(1,10) zeros(1,12) nan(1,5)]';
img = [nan(1,25) zeros(1,25) nan(1,5) zeros(1,5) nan(1,5) zeros(1,5) nan(1,10) zeros(1,5) nan(1,15) zeros(1,5) nan(1,20) zeros(1,5) nan(1,30) zeros(1,5) nan(1,50) zeros(1,60) nan(1,25)]';
img= repmat(img,1,30);
shape = img2shape(img,[0 NaN],10);
shape(:,1)=shape(:,1)+5;

fprintf('Processing stimulus.\n')
tic
s = stim_indent_shape(shape,stim_ramp(2));
fprintf('Done   '),toc

%% compute response
fprintf('Computing response.\n')
tic
ap=affpop_linear(0.2,35);           % create full affpop
%ap=affpop_linear(0.2,35,[],[],'model','GLM');           % create full affpop
ap.afferents(ap.iPC)=[];            % remove PCs
modelidx=[ap.afferents(:).idx]';    % get model indices
rc=ap.response(s);                  % compute response
fprintf('Done   '),toc


%% parse rates
allrates=rc.rate;
% loop over SA1s
for n=1:length(unique(modelidx(ap.iSA1)))
    rates_SA(:,n) = allrates(modelidx==n&ap.iSA1');
end
% loop over RAs
for n=1:length(unique(modelidx(ap.iRA)))
    rates_RA(:,n) = allrates(modelidx==n&ap.iRA');
end
%% plots
close all
stimxz=(0:350)'/10;
[~,i]=intersect(stimxz(:,1),shape(:,1));
stimxz(i,2)=1;

%%%%%% debug
slen=size(s.trace,1);
strain=zeros(slen,length(modelidx));
udyn=zeros(slen,length(modelidx));
for ii= 1:length(rc.responses)
    strain(:,ii)=rc.responses(ii).propagated_struct.stat_comp;
    udyn(:,ii)=rc.responses(ii).propagated_struct.dyn_comp;
end
%%%%% debug

f1=figure('pos',[380 560 560 420]);
subplot(3,1,[1 2]), hold on
for n=1:length(unique(modelidx(ap.iSA1)))
    idx=modelidx==n&ap.iSA1';
    plot(ap.location(idx,1),rates_SA(:,n),'b')
end
set(gca,'xticklabel',[])
ylabel('Firing rate [spk/s]')
title('SA1')

subplot(3,1,3), hold on
st=strain(round(slen/2),idx);
ud=mean(udyn(1:500,idx));

plot(stimxz(:,1),stimxz(:,2),'-','color',[0 0 0 .7])
plot(ap.location(idx,1),st/max(st),'color',[1 0 0 .7])
plot(ap.location(idx,1),ud/max(ud),'color',[0 .7 0 .7])
set(gca,'yticklabel',[])
xlabel('Position [mm]')
hleg=legend('Stim','Strain','Dyn'); 
set(hleg,'box','off','pos',get(hleg,'pos')+[.1 .1 0 0])

screenshot(f1,'figs/spatial_res_bars01_curr','pdf')
close(f1)

f2=figure('pos',[980 560 560 420]);
subplot(3,1,[1 2]), hold on
for n=1:length(unique(modelidx(ap.iRA)))
    idx=modelidx==n&ap.iRA';
    plot(ap.location(idx,1),rates_RA(:,n),'b')
end
set(gca,'xticklabel',[])
ylabel('Firing rate [spk/s]')
title('RA')

subplot(3,1,3), hold on
st=strain(round(slen/2),idx);
ud=mean(udyn(1:500,idx));

plot(stimxz(:,1),stimxz(:,2),'-','color',[0 0 0 .7])
plot(ap.location(idx,1),st/max(st),'color',[1 0 0 .7])
plot(ap.location(idx,1),ud/max(ud),'color',[0 .7 0 .7])
set(gca,'yticklabel',[])
xlabel('Position [mm]')
hleg=legend('Stim','Strain','Dyn'); 
set(hleg,'box','off','pos',get(hleg,'pos')+[.1 .1 0 0])

screenshot(f2,'figs/spatial_res_bars02_curr','pdf')
close(f2)

% screenshot(1)
% screenshot(1,'SAcomp1')
% screenshot(2,'RAcomp1')
% screenshot(0)