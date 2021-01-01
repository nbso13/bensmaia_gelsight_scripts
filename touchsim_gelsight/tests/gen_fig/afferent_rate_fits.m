function afferent_rate_fits
close all;
a = affpop_single_models();

%% Sinusoids

try
    load SINE.mat
catch
    try
        load \\bsl-somsrv1\data\processed\per\pma\L2_str\SINE.mat
    catch
        error('Unable to find SINE.mat locally or on server.')
    end
end

%      SA         RA                         PC
idx = [2 6 7 10   1 3 13 14 15 16 22 23 25   8 9 11 27];
num = length(idx);

stim_len = SINE.state_v(:,8);
idx_sine = 1:length(stim_len);

[freqs,ind_sorted] = sort(SINE.state_v(idx_sine,5),'ascend');
amps=SINE.state_v(idx_sine(ind_sorted),2);
ntypes=SINE.nType(idx);
freqs_unique = unique(freqs);
idx_sine = idx_sine(ind_sorted);
stim_len = stim_len(ind_sorted);

try
    load desiredsinepos.mat
catch
    process_sine_pos;
    load desiredsinepos.mat
end

for i=1:length(idx_sine)
    for n=1:num
        r_act(n,i) = length(SINE.ObservedSpikes{idx(n)}{idx_sine(i),1})/...
            stim_len(i);
        times_act(n,i)=SINE.ObservedSpikes{idx(n)}(idx_sine(i),1);
    end
    
    s = Stimulus(desiredsinepos{i}'/1000,[0 0],20000,0.5);
    r = a.response(s);
    
    r_sim(:,i) = r.rate;
    times_sim(:,i)={r.responses(:).spikes};
    
    %% debug
    %b=[r.responses(:).propagated_struct];
    %st=[b(:).stat_comp];
    %rmsstrain(i,:)=rms(st(:,[1 10 21]));
    
end

% compute timing specific variables
% number of spikes
n_act=cellfun(@length,times_act);
n_sim=cellfun(@length,times_sim);
nidx=[1:4 1:9 1:4];

% ISI
isi_act=cellfun(@diff,times_act,'uni',0);
isi_sim=cellfun(@diff,times_sim,'uni',0);

% median ISI in a trial
medisi_act=cellfun(@median,isi_act);
medisi_sim=cellfun(@median,isi_sim);

nn_act=~isnan(medisi_act); % non-nan indices
nn_sim=~isnan(medisi_sim);

% thresholds
for f=1:length(freqs_unique)
    idx = find(freqs==freqs_unique(f));
    medisi_act_f = medisi_act(:,idx);
    medisi_sim_f = medisi_sim(:,idx);
    for n=1:size(medisi_act,1)
        min_act_f = min(amps(idx(~isnan(medisi_act_f(n,:)))));
        if isempty(min_act_f)
            thres_act(n,f) = NaN;
        else
            thres_act(n,f) = min_act_f;
        end
        min_sim_f = min(amps(idx(~isnan(medisi_sim_f(n,:)))));
        if isempty(min_sim_f)
            thres_sim(n,f) = NaN;
        else
            thres_sim(n,f) = min_sim_f;
        end
    end
end

%% Noise

try
    load NOISE.mat
catch
    try
        load \\bsl-somsrv1\data\processed\per\pma\L2_str\NOISE.mat
    catch
        error('Unable to find NOISE.mat locally or on server.')
    end
end

%      SA        RA                         PC
idx = [2 4 5 8   1 3 11 12 13 14 18 19 20   6 7 9 21];

try
    load desirednoisepos.mat
catch
    process_noise_pos;
    load desirednoisepos.mat
end

for i=1:60
    for n=1:num
        r_act_noise(n,i) = length(NOISE.ObservedSpikes{idx(n)}{i,1});
    end
    
    s = Stimulus(reshape(desirednoisepos{i},[],1)/1000,[0 0],20000,0.5);
    r = a.response(s);
    r_resp = r.responses;
    
    r_sim_noise(:,i) = r.rate;
end

%% Diharms

try
    load DI.mat
catch
    try
        load \\bsl-somsrv1\data\processed\per\pma\L2_str\DI.mat
    catch
        error('Unable to find DI.mat locally or on server.')
    end
end

%      SA         RA                         PC
idx = [2 6 7 10   1 3 13 14 15 16 21 22 24   8 9 11 27];
num = length(idx);

try
    load desireddiharmpos.mat
catch
    process_diharm_pos;
    load desireddiharmpos.mat
end

for i=1:size(stim,1)
    stim_ind = find(DI.state_v(:,5)==stim(i,1) & DI.state_v(:,6)==stim(i,2) & DI.state_v(:,2)==stim(i,3));
    for n=1:num
        r_act_di(n,i) = length(DI.ObservedSpikes{idx(n)}{stim_ind,1})/DI.state_v(stim_ind,end);
    end
    
    s = Stimulus(desireddiharmpos{i}/1000,[0 0],20000,0.5);
    r = a.response(s);
    r_sim_di(:,i) = r.rate;
end

%% Actual vs predicted rates for sinusoids

figure(1)
set(1,'pos',[1 1 1920 1080],'PaperPosition', [.25 .25 25 25]);
for n=1:num
    subplot(4,5,n)
    hold on
    if strcmpi(a.afferents(n).class,'PC')
        F = 17;
        cc_tmp = corrcoef(r_act(n,:),r_sim(n,:),'rows','complete');
    else
        F = 8;
        cc_tmp = corrcoef(r_act(n,1:120),r_sim(n,1:120),'rows','complete');
    end
    for f=1:F
        plot(r_act(n,freqs==freqs_unique(f)),r_sim(n,freqs==freqs_unique(f)),'o','Color',[1-log10(freqs_unique(f))/3 log10(freqs_unique(f))/3 0])
    end
    unityslope;
    cc(n) = cc_tmp(2);
    title([a.afferents(n).class ': r=' num2str(cc_tmp(2))])
end

screenshot(1,'figs/afferent_rate_fits01_curr','pdf')
close(1)

%% Actual vs predicted rates for noise

figure(2)
set(2,'pos',[1 1 1920 1080],'PaperPosition', [.25 .25 25 25]);
for n=1:num
    subplot(4,5,n)
    hold on
    if strcmpi(a.afferents(n).class,'PC')
        F = 1;
    else
        F = 41;
    end
    plot(r_act_noise(n,F:end),r_sim_noise(n,F:end),'ko')
    unityslope;
    cc_tmp = corrcoef(r_act_noise(n,F:end),r_sim_noise(n,F:end),'rows','complete');
    cc(n) = cc_tmp(2);
    title([a.afferents(n).class ': r=' num2str(cc_tmp(2))])
end

screenshot(2,'figs/afferent_rate_fits02_curr','pdf')
close(2)

%% Actual vs predicted rates for diharms

figure(3)
set(3,'pos',[1 1 1920 1080],'PaperPosition', [.25 .25 25 25]);
for n=1:num
    subplot(4,5,n)
    if strcmpi(a.afferents(n).class,'PC')
        F = [1:138];
    else
        F = [1:55 76:80];
    end
    plot(r_act_di(n,F),r_sim_di(n,F),'o')
    unityslope;
    cc_tmp = corrcoef(r_act_di(n,F),r_sim_di(n,F));
    cc(n) = cc_tmp(2);
    title([a.afferents(n).class ': r=' num2str(cc_tmp(2))])
end

screenshot(3,'figs/afferent_rate_fits03_curr','pdf')
close(3)

%% Actual vs predicted thresholds for sinusoids

figure(9)
set(9,'pos',[1 1 1920 1080],'PaperPosition', [.25 .25 25 25]);
for n=1:num
    subplot(4,5,n)
    loglog(thres_act(n,:),thres_sim(n,:),'-o')
    unityslope;
    cc_tmp = corrcoef(log(thres_act(n,:)),log(thres_sim(n,:)),'rows','complete');
    cc(n) = cc_tmp(2);
    title([a.afferents(n).class ': r=' num2str(cc_tmp(2))])
end

screenshot(9,'figs/afferent_rate_fits09_curr','pdf')
close(9)

%% median ISI corr

figure(4)
v=logspace(-3,0,200);
clear ax
for ii = 1:18
    ax(ii) = subplot(4,5,ii);
end
set(ax,'nextplot','add')
for ii=1:num
    act=medisi_act(ii,:)';    sim=medisi_sim(ii,:)';
    badnan=isnan(act)|isnan(sim);
    ab=[act(~badnan) ones(size(act(~badnan)))]\sim(~badnan);
    scatter(ax(ii),act,sim,[],freqs(:,1));
    plot(ax(ii),v,v,v,ab(1)*v+ab(2));
    title(ax(ii),[ntypes{ii} ' ' num2str(nidx(ii)) ':   '...
        num2str(corr(act(~badnan),sim(~badnan)))])
end
set(ax(18),'clim',[1 1e3],'visible','off'); axes(ax(18));colorbar;
xlabel(ax(17),'Actual ISI [s]');
ylabel(ax(11),'Simulated ISI [s]');
set(ax,'xscale','log','yscale','log','xlim',[1e-3 1e0],'ylim',[1e-3 1e0])
set(gcf,'pos',[1 1 1920 1080]);

screenshot(gcf,'figs/afferent_rate_fits04_curr','pdf')
close(4)


%% ISI 2D plot (freq/amp)
clear ax axh

nneur=[4 5 4 4];
nnumb=[ones(1,4) ones(1,5)*2 ones(1,4)*3 ones(1,4)*4];
nidxs=[1:4 1:5 1:4 1:4];

% plot all
for ff=1:4
    figure(ff)
    set(ff,'pos',[1    61   400*nneur(ff)   750])
    for ii = 1:(2*(nneur(ff)+1))
        axh{ff}(ii) = subplot(2,nneur(ff)+1,ii);
    end
    for ii=1:length(axh{ff})
        p=get(axh{ff}(ii),'pos');
        p(1)=(p(1)-.5)*1.20+.5;
        set(axh{ff}(ii),'pos',p);
    end
    for nn=1:nneur(ff)
        ax{ff,nn}=axh{ff}([nn nn+nneur(ff)+1]);
    end
    colormap(ff,flipud(jet))
end
tri=delaunay(freqs,amps);

for id=1:num
    f=nnumb(id);
    n=nidxs(id);
    set(ax{f,n},'nextplot','add');
    % plot empty circles
    scatter(ax{f,n}(1),freqs,amps,40,[.7 .7 .7])
    scatter(ax{f,n}(2),freqs,amps,40,[.7 .7 .7])
    % plot coloredfilled dots
    idx1=find(~isnan(medisi_act(id,:))& medisi_act(id,:)<1);
    scatter(ax{f,n}(1),freqs(idx1),amps(idx1),40,log10(medisi_act(id,idx1)),'filled')
    idx2=find(~isnan(medisi_sim(id,:)));
    scatter(ax{f,n}(2),freqs(idx2),amps(idx2),40,log10(medisi_sim(id,idx2)),'filled')
    % draw a line for thresholds
    plot(ax{f,n}(1),freqs_unique,thres_act(id,:),'k')
    plot(ax{f,n}(2),freqs_unique,thres_sim(id,:),'k')
    
    
    set(ax{f,n},'xscale','log','yscale','log','xlim',[1 1000],...
        'ylim',[.0625 500],'clim',[-3 0],...
        'xtick',[1 10 100 1000],'ytick',[1 10 100],'xminorgrid','off',...
        'yminorgrid','off');
    xlabel(ax{f,n}(1),'Frequency [Hz]')
    ylabel(ax{f,n}(1),'Amplitude [um]')
    title(ax{f,n}(1),['Median ISI (' ntypes{id} num2str(nidx(id)) ') : ACTUAL']);
    xlabel(ax{f,n}(2),'Frequency [Hz]')
    ylabel(ax{f,n}(2),'Amplitude [um]')
    title(ax{f,n}(2),'SIMULATED');
end

for ff=1:4
    h(1)=colorbar(axh{ff}(end/2));
    h(2)=colorbar(axh{ff}(end));
    set(axh{ff}([end/2 end]),'visible','off','clim',[-3 0])
    set(h,'ticks',[-3 -2 -1 0],'ticklabels',[.001 .01 .1 1]);
    
    h(1).Label.String='ISI [s]';
    h(2).Label.String='ISI [s]';
    
    screenshot(ff,['figs/afferent_rate_fits0' num2str(ff+4) '_curr'],'pdf')
    close(ff)
end

%% Rate-intensity functions
clear ax axh

nneur=[4 5 4 4];
nnumb=[ones(1,4) ones(1,5)*2 ones(1,4)*3 ones(1,4)*4];
nidxs=[1:4 1:5 1:4 1:4];

% plot all
for ff=1:4
    figure(ff)
    set(ff,'pos',[1    61   400*nneur(ff)   750])
    for ii = 1:(2*(nneur(ff)+1))
        axh{ff}(ii) = subplot(2,nneur(ff)+1,ii);
    end
    for ii=1:length(axh{ff})
        p=get(axh{ff}(ii),'pos');
        p(1)=(p(1)-.5)*1.20+.5;
        set(axh{ff}(ii),'pos',p);
    end
    for nn=1:nneur(ff)
        ax{ff,nn}=axh{ff}([nn nn+nneur(ff)+1]);
    end
end

for id=1:num
    f=nnumb(id);
    n=nidxs(id);
    set(ax{f,n},'nextplot','add');
    % plot RI functions for all frequencies
    for fff=1:length(freqs_unique)
        idx = freqs==freqs_unique(fff);
        plot(ax{f,n}(1),amps(idx),r_act(id,idx))
        plot(ax{f,n}(2),amps(idx),r_sim(id,idx))
    end
    
    set(ax{f,n},'xscale','log','xlim',[1 1000],...
        'xtick',[1 10 100 1000],'xminorgrid','off',...
        'yminorgrid','off');
    xlabel(ax{f,n}(1),'Frequency [Hz]')
    ylabel(ax{f,n}(1),'Firing rate [Hz]')
    title(ax{f,n}(1),['Rate-intensity (' ntypes{id} num2str(nidx(id)) ') : ACTUAL']);
    xlabel(ax{f,n}(2),'Frequency [Hz]')
    ylabel(ax{f,n}(2),'Firing rate [Hz]')
    title(ax{f,n}(2),'SIMULATED');
end

for ff=1:4
    screenshot(ff,['figs/afferent_rate_fits' num2str(ff+9) '_curr'],'pdf')
    close(ff)
end


