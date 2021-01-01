% spike_timing.m

addpath helper_functions/

cols=affcol;

%% RA entrainment

ap = affpop_single_models([],[],'noisy',false);
ap.afferents(~ap.iRA)=[]; % only RA's
naff=length(ap.afferents);

frq=40;
rad=1;
pre_indent=.55;
ampRA = linspace(0,0.3,50);
namp=length(ampRA);

resp_rateRA=zeros(namp,naff);
phase_angleRA=NaN*zeros(namp,naff);
phase_angle2RA=NaN*zeros(namp,naff);
v_strengthRA=NaN*zeros(namp,naff);
for ii=1:namp
    stim=stim_sine(frq,ampRA(ii),[],[],[],[],[],rad,pre_indent);
    rc(ii)=ap.response(stim);
    resp_rateRA(ii,:)=rc(ii).rate;

    for jj=1:naff
        sp_tmp = rc(ii).responses(jj).spikes;
        if ~isempty(sp_tmp)
            
            hh = histc(sp_tmp,0:0.025:1);
            hh_idx = find(hh>1);
            del = [];
            idx2 = [];
            for kk=1:length(hh_idx)
                del(end+1:end+hh(hh_idx(kk))-1) = sum(hh(1:hh_idx(kk)-1))+2:sum(hh(1:hh_idx(kk)));
                idx2(end+1) = sum(hh(1:hh_idx(kk)-1))+2;
            end
            sp2_tmp = sp_tmp(idx2);
            sp_tmp(del) = [];
            phase_angleRA(ii,jj) = circ_mean(mod(sp_tmp,0.025)*2*pi/0.025,[],2);
            v_strengthRA(ii,jj) = vector_strength(mod(sp_tmp,0.025)*2*pi/0.025);
            
            if ~isempty(sp2_tmp)
                phase_angle2RA(ii,jj) = mean(mod(sp2_tmp,0.025)*360/0.025);
            end
        end
    end
end
phase_angleRA = unwrap(phase_angleRA);

for i=1:9
    idx = resp_rateRA(:,i)>35 & resp_rateRA(:,i)<45;
    I1idxRA(i) = find(idx,1);
    I1RA(i) = ampRA(I1idxRA(i));
end

amp_rat = logspace(-0.3,0.7,25);
for i=1:9
    paRA(:,i) = interp1(ampRA(:)./I1RA(i),phase_angleRA(:,i)-mean(phase_angleRA(end-10:end,i)),amp_rat);
    vs_thresRA(i) = v_strengthRA(I1idxRA(i),i);
end

%% PC entrainment

ap = affpop_single_models();
ap.afferents(~ap.iPC)=[]; % only PC's
naff=length(ap.afferents);

frq=200;
rad=1;
pre_indent=.55;
ampPC = linspace(0,0.025,50);
namp=length(ampPC);

resp_ratePC = zeros(namp,naff);
phase_anglePC = NaN*zeros(namp,naff);
phase_angle2PC = NaN*zeros(namp,naff);
V_strengthPC = NaN*zeros(namp,naff);
for ii=1:namp
    stim=stim_sine(frq,ampPC(ii),[],[],[],[],[],rad,pre_indent);
    rcPC(ii)=ap.response(stim);
    resp_ratePC(ii,:)=rcPC(ii).rate;

    for jj=1:naff
        sp_tmp = rcPC(ii).responses(jj).spikes;
        if ~isempty(sp_tmp)
            
            hh = histc(sp_tmp,0:0.025:1);
            hh_idx = find(hh>1);
            del = [];
            idx2 = [];
            for kk=1:length(hh_idx)
                del(end+1:end+hh(hh_idx(kk))-1) = sum(hh(1:hh_idx(kk)-1))+2:sum(hh(1:hh_idx(kk)));
                idx2(end+1) = sum(hh(1:hh_idx(kk)-1))+2;
            end
            sp2_tmp = sp_tmp(idx2);
            sp_tmp(del) = [];
            phase_anglePC(ii,jj) = circ_mean(mod(sp_tmp,0.025)*2*pi/0.025,[],2);
            v_strengthPC(ii,jj) = vector_strength(mod(sp_tmp,0.025)*2*pi/0.025);
            
            if ~isempty(sp2_tmp)
                phase_angle2PC(ii,jj) = mean(mod(sp2_tmp,0.025)*360/0.025);
            end
        end
    end
end
phase_anglePC = unwrap(phase_anglePC);

for i=1:4
    idx = resp_ratePC(:,i)>180 & resp_ratePC(:,i)<220;
    I1idxPC(i) = find(idx,1);
    I1PC(i) = ampPC(I1idxPC(i));
end

amp_rat = logspace(-0.3,1,25);
for i=1:4
    paPC(:,i) = interp1(ampPC(:)./I1PC(i),phase_anglePC(:,i)-mean(phase_anglePC(end-10:end,i)),amp_rat);
    vs_thresPC(i) = v_strengthPC(I1idxPC(i),i);
end

%%

figure(1)
set(gcf,'pos',[0 0 240 474]);

subplot(211)
semilogx(amp_rat,nanmean(paRA,2)/2/pi*360,'Color',cols(2,:),'LineWidth',1.5)
hold on
semilogx(amp_rat,nanmean(paPC,2)/2/pi*360,'Color',cols(3,:),'LineWidth',1.5)
line([1 1],[-25 150],'Color','k')
ylim([-10 150])
xlim([.5 10])
box off
ylabel('Phase angle [°]')
xlabel('Amplitude/tuning point')

subplot(212)
plot(1+randn(1,9)/7.5,vs_thresRA,'o','Color',cols(2,:))
hold on
line([0.75 1.25],[mean(vs_thresRA) mean(vs_thresRA)],'LineWidth',1.5,'Color','k')
plot(2+randn(1,4)/7.5,vs_thresPC,'o','Color',cols(3,:))
hold on
line([1.75 2.25],[mean(vs_thresPC) mean(vs_thresPC)],'LineWidth',1.5,'Color','k')
ylim([0 1])
set(gca,'xtick',[1 2],'xticklabel',{'RA','PC'})
box off
ylabel('Vector strength')

screenshot(1,'spike_timing1','pdf')
close(1)

%%

ex_id = 6;

X1 = csvread('spike_timing2a.csv',1);
X2 = csvread('spike_timing2b.csv',1);

figure(2)
set(gcf,'pos',[0 0 240 474]);

subplot(411)
plot(X2(:,1)/1000,X2(:,2),'k','LineWidth',2)
xlim([0 0.3])
box off
ylabel('Phase [°]')

subplot(412)
plot(X1(:,1)/1000,X1(:,2),'k','LineWidth',2)
xlim([0 0.3])
box off
ylabel('Impulse/cycle')
xlabel('Indentation depth [mm]')

subplot(413)
plot(ampRA,mod(phase_angleRA(:,ex_id)-mean(phase_angleRA(end-10:end,ex_id))+0.85,2*pi)/2/pi*360,'LineWidth',2)
box off
ylabel('Phase [°]')
xlim([0 .2])
ylim([0 150])

subplot(414)
plot(ampRA,resp_rateRA(:,ex_id)/40,'LineWidth',2)
box off
ylim([0 2])
xlim([0 .2])
ylabel('Impulse/cycle')
xlabel('Indentation depth [mm]')

screenshot(2,'spike_timing2','pdf')
close(2)
