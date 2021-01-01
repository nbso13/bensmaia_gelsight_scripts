function decay_dist_RA
% clear; close all; clc
%% create affpop
locs = logspace(-1,1,10);
nlocs=length(locs);

ap = AfferentPopulation;
for l=1:nlocs
    a_sub = affpop_single_models([locs(l) 0]);
    a_sub.afferents(~a_sub.iRA)=[]; % only RA's
    ap.afferents = [ap.afferents a_sub.afferents];
end
naff=length(ap.afferents);

%% create stimuli and compute resp

frq=40;
rad=1;
pre_indent=.55;
%pre_indent=0;
amp = linspace(0,1,60);
namp=length(amp);

resp_rate=zeros(namp,naff);
% resp_proprms=zeros(namp,naff); % FOR DEBUG
phase_angle=NaN*zeros(namp,naff);
phase_angle2=NaN*zeros(namp,naff);
for ii=1:namp
    stim=stim_sine(frq,amp(ii),[],[],[],[],[],rad,pre_indent);
    rc(ii)=ap.response(stim);
    resp_rate(ii,:)=rc(ii).rate;
%     for jj=1:naff
%         resp_proprms(ii,jj)=rms(rc(ii).responses(jj).propagated_struct.stat_comp);
%     end

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
            phase_angle(ii,jj) = mean(mod(sp_tmp,0.025)*360/0.025);
            if ~isempty(sp2_tmp)
                phase_angle2(ii,jj) = mean(mod(sp2_tmp,0.025)*360/0.025);
            end
        end
    end
    
    fprintf('.')
end
fprintf('\n')

%%
cols=affcol;

% tic
idx_phaseloc=zeros(naff,1);
idx_resp=zeros(naff,1);
for ii=1:naff
    % find minimal (first) ampl showing phase locked resp (.9 spikes/cycle)
    truc=find(resp_rate(:,ii)>frq*.9,1);
    if(isempty(truc)),       truc=namp;   end
    idx_phaseloc(ii)=truc;
    
    % find minimal (first) ampl showing a response (.1 spikes/s)
    truc=find(resp_rate(:,ii)>.1,1);
    if(isempty(truc)),       truc=namp;   end
    idx_resp(ii)=truc;
end
% toc

% reshape phaseloc resp to matrix nmodels x nlocs
amp_phaseloc=reshape(amp(idx_phaseloc),naff/nlocs,nlocs);
amp_resp=reshape(amp(idx_resp),naff/nlocs,nlocs);


%% basic plots (for debug)
debug=0;
if(debug)
    % rates for each neuron model
    ax=subplot_ax(4,4);
    for ii=1:14
        idx=(ii-1)*10+1:ii*10;
        plot(ax(ii),locs,resp_rate(:,idx));
    end
    legend(ax(14),cellstr(num2str(amp')))
end
%% Fig 7 Johnson 1973

close all
figure(1)
set(1,'defaultlinelinewidth',2)
set(1,'DefaultAxesFontSize',14)
errorbar(locs,nanmean(amp_phaseloc),nansem(amp_phaseloc,0,1),'color',cols(2,:))
set(gca,'xscale','log','yscale','log','xlim',[.1 10],'ylim',[.01 1])
xlabel('Distance from stimulus [mm]')
ylabel('Minimal ampl. eliciting phase loc response [mm]')
title('Simulation Fig 7 (Johnson 1974)','fontsize',16)

screenshot(1,'figs/decay_dist_RA01_curr','pdf')
close(1)

%% Fig 6 Johnson 1973
figure(2)
set(2,'defaultlinelinewidth',2)
set(2,'DefaultAxesFontSize',14)
plot(locs,amp_resp,'.','markersize',25,'color',cols(2,:))
set(gca,'xscale','log','yscale','log','ylim',[.001 1],'xlim',[1 14])
xlabel('Distance from stimulus [mm]')
ylabel('Minimal ampl. eleciting response [mm]')
title('Simulation Fig 6 (Johnson 1974)','fontsize',16)

screenshot(2,'figs/decay_dist_RA02_curr','pdf')
close(2)

%% Fig 3 Johnson 1973

figure(3)
for i=1:9
    subplot(4,4,i)
    hold on
    plot(amp(1:20),resp_rate(1:20,i)/40,'k')
    idx = resp_rate(1:20,i)>35 & resp_rate(1:20,i)<45;
    line(amp([find(idx,1) find(idx,1)]),[0 3],'Color','k')
    line(amp([find(idx,1,'last') find(idx,1,'last')]),[0 3],'Color','k')
    plot(amp(1:20),phase_angle(1:20,i)/360*3)
    plot(amp(1:20),phase_angle2(1:20,i)/360*3)
    xlim([0 0.3])
    ylim([0 3])
    title(['RA ' num2str(i)])
    set(gca,'ygrid','on')
end

screenshot(3,'figs/decay_dist_RA03_curr','pdf')
close(3)
