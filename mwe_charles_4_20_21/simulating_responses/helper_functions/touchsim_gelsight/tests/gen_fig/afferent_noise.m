function afferent_noise
%clear all;close all;clc
a = affpop_single_models([],[],'noisy',true);
tc = [repmat(250,1,17) repmat(100,1,4)];

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
        times_act(n,i) = {SINE.ObservedSpikes{idx(n)}(idx_sine(i),1:5)};
    end
    
    s = Stimulus(desiredsinepos{i}'/1000,[0 0],20000,0.5);
    for rep=1:5
        r_resp{rep} = a.response(s).responses;
    end
    
    for n=1:num
        times_sim(n,i) = {{r_resp{1}(n).spikes r_resp{2}(n).spikes r_resp{3}(n).spikes r_resp{4}(n).spikes r_resp{5}(n).spikes}};
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
        times_act_noise(n,i) = {NOISE.ObservedSpikes{idx(n)}(i,1:5)};
    end
    
    s = Stimulus(reshape(desirednoisepos{i},[],1)/1000,[0 0],20000,0.5);
    for rep=1:5
        r_resp{rep} = a.response(s).responses;
    end
    
    for n=1:num
        times_sim_noise(n,i) = {{r_resp{1}(n).spikes r_resp{2}(n).spikes r_resp{3}(n).spikes r_resp{4}(n).spikes r_resp{5}(n).spikes}};
    end
end

%% calculate noise level

for n=1:num
    for i=1:length(idx_sine)
        c = 1;
        for r1=1:5
            for r2=r1+1:5
                edges = 0:1/5000:length(desiredsinepos{i})/20000;
                sbin1 = histc(times_act{n,i}{r1}-0.05,edges)';
                sbin2 = histc(times_act{n,i}{r2}-0.05,edges)';
                dists_act_tmp(c) = VRdist_pairwise(sbin1(:),sbin2(:),tc(n));
                
                sbin1 = histc(times_sim{n,i}{r1},edges)';
                sbin2 = histc(times_sim{n,i}{r2},edges)';
                dists_sim_tmp(c) = VRdist_pairwise(sbin1(:),sbin2(:),tc(n));
                
                c = c + 1;
            end
        end
        dists_act(n,i) = nanmean(dists_act_tmp);
        dists_sim(n,i) = nanmean(dists_sim_tmp);
    end
    
    for i=1:60
        c = 1;
        for r1=1:5
            for r2=r1+1:5
                edges = 0:1/5000:length(desirednoisepos{i})/20000;
                sbin1 = histc(times_act_noise{n,i}{r1}-0.5,edges)';
                sbin2 = histc(times_act_noise{n,i}{r2}-0.5,edges)';
                dists_act_tmp(c) = VRdist_pairwise(sbin1(:),sbin2(:),tc(n));
                
                sbin1 = histc(times_sim_noise{n,i}{r1},edges)';
                sbin2 = histc(times_sim_noise{n,i}{r2},edges)';
                dists_sim_tmp(c) = VRdist_pairwise(sbin1(:),sbin2(:),tc(n));
                
                c = c + 1;
            end
        end
        dists_act_noise(n,i) = nanmean(dists_act_tmp);
        dists_sim_noise(n,i) = nanmean(dists_sim_tmp);
    end
end

%% Noise plots for sinusoids

figure(5)
set(5,'pos',[1 1 1920 1080],'PaperPosition', [.25 .25 25 25]);
for n=1:num
    subplot(4,5,n)
    hold on
    cc_tmp = corrcoef(dists_act(n,:),dists_sim(n,:),'rows','complete');
    
    plot(dists_act(n,:),dists_sim(n,:),'ko')

    unityslope;
    cc(n) = cc_tmp(2);
    title([a.afferents(n).class ': r=' num2str(cc_tmp(2))])
end

screenshot(5,'figs/afferent_noise05_curr','pdf')
close(5)

%% Noise plots for noise

figure(6)
set(6,'pos',[1 1 1920 1080],'PaperPosition', [.25 .25 25 25]);
for n=1:num
    subplot(4,5,n)
    hold on
    cc_tmp = corrcoef(dists_act_noise(n,:),dists_sim_noise(n,:),'rows','complete');
    
    plot(dists_act_noise(n,:),dists_sim_noise(n,:),'ko')

    unityslope;
    cc(n) = cc_tmp(2);
    title([a.afferents(n).class ': r=' num2str(cc_tmp(2))])
end

screenshot(6,'figs/afferent_noise06_curr','pdf')
close(6)

%% Spike train plots

clear ax axh

nneur=[4 5 4 4];
nnumb=[ones(1,4) ones(1,5)*2 ones(1,4)*3 ones(1,4)*4];
nidxs=[1:4 1:5 1:4 1:4];
nidx=[1:4 1:9 1:4];

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
    
    switch ntypes{id}
        case 'SA'
            idx_noise = 58;
            idx_sine = 90;
            dur_noise = 1;
            dur_sine = .1;
        case 'RA'
            idx_noise = 57;
            idx_sine = 86;
            dur_sine = .1;
        case 'PC'
            idx_noise = 16;
            idx_sine = 127;
            dur_noise = .2;
            dur_sine = .1;
    end
    
    % plot spike trains for noise stimulus
    plot_spikes(times_act_noise{id,idx_noise},'col','k','par',ax{f,n}(1),'neuron_offset',5,'time_o',-.5)
    plot_spikes(times_sim_noise{id,idx_noise},'col','b','par',ax{f,n}(1))

    % plot spike trains for sine stimulus
    plot_spikes(times_act{id,idx_sine},'col','k','par',ax{f,n}(2),'neuron_offset',5,'time_o',-.05)
    plot_spikes(times_sim{id,idx_sine},'col','b','par',ax{f,n}(2))
    
    set(ax{f,n}(1),'xlim',[0 dur_noise])
    set(ax{f,n}(2),'xlim',[0 dur_sine])
    xlabel(ax{f,n}(2),'Time (s)')
    title(ax{f,n}(1),['Spike trains (' ntypes{id} num2str(nidx(id)) ') : Noise']);
    title(ax{f,n}(2),'Sine');
end

for ff=1:4
    screenshot(ff,['figs/afferent_noise0' num2str(ff) '_curr'],'pdf')
    close(ff)
end
