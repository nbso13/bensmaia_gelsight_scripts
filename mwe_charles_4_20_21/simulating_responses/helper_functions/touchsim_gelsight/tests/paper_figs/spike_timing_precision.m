% spike_timing_precision.m
shouldUseGLM = 1;

RAfreq = 40;
PCfreq = 300;

tuning_thresh = 0.9;

addpath helper_functions/

a = affpop_single_models([],[],'noisy',true);
if(exist('shouldUseGLM','var') && shouldUseGLM)
    a = affpop_single_models([],[],'Model','GLM','noisy',true);
end
tc = [repmat(250,1,17) repmat(100,1,4)];

%% Sinusoids
load SINE.mat

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

n_act = nan(num,length(idx_sine),5);
n_sim = nan(num,length(idx_sine),5);
vs_act = nan(num,length(idx_sine),5);
vs_sim = nan(num,length(idx_sine),5);

for i=1:length(idx_sine)
    for n=1:num
        times_act(n,i) = {SINE.ObservedSpikes{idx(n)}(idx_sine(i),1:5)};
        
        for rep=1:5
            if length(times_act{n,i}{rep})>=freqs(i)/10
                sp_tmp = times_act{n,i}{rep};
                hh = histc(sp_tmp,0:(1/freqs(i)):1);
                hh_idx = find(hh>1);
                del = [];
                for kk=1:length(hh_idx)
                    del(end+1:end+hh(hh_idx(kk))-1) = sum(hh(1:hh_idx(kk)-1))+2:sum(hh(1:hh_idx(kk)));
                end
                sp_tmp(del) = [];
                n_act(n,i,rep) = length(sp_tmp)/freqs(i)/stim_len(i);
                vs_act(n,i,rep) = vector_strength(mod(sp_tmp,1/freqs(i))*2*pi/(1/freqs(i)));
                
            end
        end
        
    end
    
    s = Stimulus(desiredsinepos{i}'/1000,[0 0],20000,0.5);
    for rep=1:5
        r_resp{rep} = a.response(s).responses;
    end
    
    for n=1:num
        times_sim(n,i) = {{r_resp{1}(n).spikes r_resp{2}(n).spikes r_resp{3}(n).spikes r_resp{4}(n).spikes r_resp{5}(n).spikes}};
        
        for rep=1:5
            if length(times_sim{n,i}{rep})>=freqs(i)/10
                sp_tmp = times_sim{n,i}{rep};
                hh = histc(sp_tmp,0:(1/freqs(i)):1);
                hh_idx = find(hh>1);
                del = [];
                for kk=1:length(hh_idx)
                    del(end+1:end+hh(hh_idx(kk))-1) = sum(hh(1:hh_idx(kk)-1))+2:sum(hh(1:hh_idx(kk)));
                end
                sp_tmp(del) = [];
                n_sim(n,i,rep) = length(sp_tmp)/freqs(i)/stim_len(i);
                vs_sim(n,i,rep) = vector_strength(mod(sp_tmp,1/freqs(i))*2*pi/(1/freqs(i)));
                
            end
        end
        
    end
end

% RA vector strength at tuning point
idx = find(freqs==RAfreq);
vsm_act = mean(vs_act(a.iRA,idx,:),3);
nm_act = mean(n_act(a.iRA,idx,:),3);
vsm_sim = mean(vs_sim(a.iRA,idx,:),3);
nm_sim = mean(n_sim(a.iRA,idx,:),3);

vsRA_act=nan(size(vsm_act,1),1);
vsRA_sim=nan(size(vsm_sim,1),1);
for i=1:size(vsm_act,1)
    vsRA_act(i) = mean(vsm_act(i,find(nm_act(i,:)>=tuning_thresh)));
    vsRA_sim(i) = mean(vsm_sim(i,find(nm_sim(i,:)>=tuning_thresh)));
end

% PC vector strength at tuning point
idx = find(freqs==PCfreq);
vsm_act = mean(vs_act(a.iPC,idx,:),3);
nm_act = mean(n_act(a.iPC,idx,:),3);
vsm_sim = mean(vs_sim(a.iPC,idx,:),3);
nm_sim = mean(n_sim(a.iPC,idx,:),3);


vsPC_act=nan(size(vsm_act,1),1);
vsPC_sim=nan(size(vsm_sim,1),1);
for i=1:size(vsm_act,1)
    vsPC_act(i) = mean(vsm_act(i,find(nm_act(i,:)>=tuning_thresh)));
    vsPC_sim(i) = mean(vsm_sim(i,find(nm_sim(i,:)>=tuning_thresh)));
end

%% Noise

load NOISE.mat

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

%% Spike train plots

%ex_id = [1 5 14];
ex_id = [4 12 17];

nneur = 3;
nnumb = ones(1,3);
nidxs = 1:3;
nidx= [1 1 1];

% plot all
figure(1)
set(1,'pos',[0 0 600 250])
for ii = 1:6
    axh{1}(ii) = subplot(2,3,ii);
end
for ii=1:length(axh{1})
    p=get(axh{1}(ii),'pos');
    p(1)=(p(1)-.5)*1.20+.5;
    set(axh{1}(ii),'pos',p);
end
for nn=1:nneur(1)
    ax{1,nn}=axh{1}([nn nn+nneur(1)]);
end

n = 1;
f = 1;
for id=ex_id
    set(ax{f,n},'nextplot','add');
    
    switch ntypes{id}
        case 'SA'
            idx_noise = 58;
            idx_sine = 90;
            dur_noise = 1;
            dur_sine = .1;
            col = affcol(1);
        case 'RA'
            idx_noise = 57;
            idx_sine = 86;
            dur_noise = 1;
            dur_sine = .1;
            col = affcol(2);
        case 'PC'
            idx_noise = 16;
            idx_sine = 127;
            dur_noise = .15;
            dur_sine = .1;
            col = affcol(3);
    end
    
    % plot spike trains for noise stimulus
    plot_spikes(times_act_noise{id,idx_noise},'col','k','par',ax{f,n}(1),'neuron_offset',5,'time_o',-.5)
    plot_spikes(times_sim_noise{id,idx_noise},'col',col,'par',ax{f,n}(1))

    % plot spike trains for sine stimulus
    plot_spikes(times_act{id,idx_sine},'col','k','par',ax{f,n}(2),'neuron_offset',5,'time_o',-.05)
    plot_spikes(times_sim{id,idx_sine},'col',col,'par',ax{f,n}(2))
    
    set(ax{f,n}(1),'xlim',[0 dur_noise])
    set(ax{f,n}(2),'xlim',[0 dur_sine])
    xlabel(ax{f,n}(2),'Time (s)')
    title(ax{f,n}(1),['Spike trains (' ntypes{id} ') : Noise']);
    title(ax{f,n}(2),'Sine');
    
    n = n+1;
end

if(~shouldUseGLM)
    screenshot(1,'spike_timing_precision1','pdf')
else
    screenshot(1,'spike_timing_precision1_GLM','pdf')
end
if(screenshot)
    close(1)
end

%% plot vector strengths

figure(2)
set(gcf,'pos',[0 0 240 474]);

subplot(211)

subplot(212)
hold on
plot(randn(1,9)/7.5,vsRA_act,'o','Color','k')
line([-.3 .3],[mean(vsRA_act) mean(vsRA_act)],'LineWidth',1.5,'Color','k')
plot(1+randn(1,9)/7.5,vsRA_sim,'o','Color',affcol(2))
line([0.7 1.3],[mean(vsRA_sim) mean(vsRA_sim)],'LineWidth',1.5,'Color','k')

plot(2.5+randn(1,4)/7.5,vsPC_act,'o','Color','k')
line([2.2 2.8],[mean(vsPC_act) mean(vsPC_act)],'LineWidth',1.5,'Color','k')
plot(3.5+randn(1,4)/7.5,vsPC_sim,'o','Color',affcol(3))
line([3.2 3.8],[mean(vsPC_sim) mean(vsPC_sim)],'LineWidth',1.5,'Color','k')
ylim([0 1])
xlim([-.4 3.9])
set(gca,'xtick',[.5 3],'xticklabel',{'RA','PC'})
box off
ylabel('Vector strength')

if(~shouldUseGLM)
    screenshot(2,'spike_timing_precision2','pdf')
else
    screenshot(1,'spike_timing_precision2_GLM','pdf')
end
if(screenshot)
close(2)
end