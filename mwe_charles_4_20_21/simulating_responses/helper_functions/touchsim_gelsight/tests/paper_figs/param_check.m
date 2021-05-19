% param_check.m

cols=affcol;

pp1_8 = load('ppPC1_8_0.mat');
pp2_8 = load('ppPC2_8_0.mat');
pp3_8 = load('ppPC3_8_0.mat');
pp4_8 = load('ppPC4_8_0.mat');

pp1_10 = load('ppPC1_10_25.mat');
pp2_10 = load('ppPC2_10_25.mat');
pp3_10 = load('ppPC3_10_25.mat');
pp4_10 = load('ppPC4_10_25.mat');

pp1_1112 = load('ppPC1_1112_0.mat');
pp2_1112 = load('ppPC2_1112_0.mat');
pp3_1112 = load('ppPC3_1112_0.mat');
pp4_1112 = load('ppPC4_1112_0.mat');


a_opt = AfferentPopulation();
a_opt.afferents = [Afferent('PC','idx',1) Afferent('PC','idx',2) Afferent('PC','idx',3) Afferent('PC','idx',4)];

a_test8 = AfferentPopulation();
a_test8.afferents = [Afferent('PC','parameters',pp1_8.pp) Afferent('PC','parameters',pp2_8.pp) Afferent('PC','parameters',pp3_8.pp) Afferent('PC','parameters',pp4_8.pp)];

a_test10 = AfferentPopulation();
a_test10.afferents = [Afferent('PC','parameters',pp1_10.pp) Afferent('PC','parameters',pp2_10.pp) Afferent('PC','parameters',pp3_10.pp) Afferent('PC','parameters',pp4_10.pp)];

a_test1112 = AfferentPopulation();
a_test1112.afferents = [Afferent('PC','parameters',pp1_1112.pp) Afferent('PC','parameters',pp2_1112.pp) Afferent('PC','parameters',pp3_1112.pp) Afferent('PC','parameters',pp4_1112.pp)];


%%

load SINE.mat

% PC afferents
idx = [8 9 11 27];
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
    
    rr_opt = a_opt.response(s);
    r_opt(:,i) = rr_opt.rate;
    times_opt(:,i)={rr_opt.responses(:).spikes};
    
    rr_test8 = a_test8.response(s);
    r_test8(:,i) = rr_test8.rate;
    times_test8(:,i)={rr_test8.responses(:).spikes};
    
    rr_test10 = a_test10.response(s);
    r_test10(:,i) = rr_test10.rate;
    times_test10(:,i)={rr_test10.responses(:).spikes};
    
    rr_test1112 = a_test1112.response(s);
    r_test1112(:,i) = rr_test1112.rate;
    times_test1112(:,i)={rr_test1112.responses(:).spikes};
    
end

% compute timing specific variables
% number of spikes
n_act=cellfun(@length,times_act);
n_opt=cellfun(@length,times_opt);
n_test8=cellfun(@length,times_test8);
n_test10=cellfun(@length,times_test10);
n_test1112=cellfun(@length,times_test1112);
nidx=[1:4 1:9 1:4];

% ISI
isi_act=cellfun(@diff,times_act,'uni',0);
isi_opt=cellfun(@diff,times_opt,'uni',0);
isi_test8=cellfun(@diff,times_test8,'uni',0);
isi_test10=cellfun(@diff,times_test10,'uni',0);
isi_test1112=cellfun(@diff,times_test1112,'uni',0);

% median ISI in a trial
medisi_act=cellfun(@median,isi_act);
medisi_opt=cellfun(@median,isi_opt);
medisi_test8=cellfun(@median,isi_test8);
medisi_test10=cellfun(@median,isi_test10);
medisi_test1112=cellfun(@median,isi_test1112);

nn_act=~isnan(medisi_act); % non-nan indices
nn_opt=~isnan(medisi_opt);
nn_test8=~isnan(medisi_test8);
nn_test10=~isnan(medisi_test10);
nn_test1112=~isnan(medisi_test1112);

% thresholds
for f=1:length(freqs_unique)
    idx = find(freqs==freqs_unique(f));
    medisi_act_f = medisi_act(:,idx);
    medisi_opt_f = medisi_opt(:,idx);
    medisi_test8_f = medisi_test8(:,idx);
    medisi_test10_f = medisi_test10(:,idx);
    medisi_test1112_f = medisi_test1112(:,idx);
    
    for n=1:size(medisi_act,1)
        min_act_f = min(amps(idx(~isnan(medisi_act_f(n,:)))));
        if isempty(min_act_f)
            thres_act(n,f) = NaN;
        else
            thres_act(n,f) = min_act_f;
        end
        
        min_opt_f = min(amps(idx(~isnan(medisi_opt_f(n,:)))));
        if isempty(min_opt_f)
            thres_opt(n,f) = NaN;
        else
            thres_opt(n,f) = min_opt_f;
        end
        
        min_test8_f = min(amps(idx(~isnan(medisi_test8_f(n,:)))));
        if isempty(min_test8_f)
            thres_test8(n,f) = NaN;
        else
            thres_test8(n,f) = min_test8_f;
        end
        
        min_test10_f = min(amps(idx(~isnan(medisi_test10_f(n,:)))));
        if isempty(min_test10_f)
            thres_test10(n,f) = NaN;
        else
            thres_test10(n,f) = min_test10_f;
        end
        
        min_test1112_f = min(amps(idx(~isnan(medisi_test1112_f(n,:)))));
        if isempty(min_test1112_f)
            thres_test1112(n,f) = NaN;
        else
            thres_test1112(n,f) = min_test1112_f;
        end
        
    end
end

%% generate summary figure


figure(1)
set(1,'pos',[1 1 700 500],'PaperPosition', [.25 .25 25 25]);

% no post-spike inhibition
n_ex = 3;
f_ex = [200 250 300];
for f=1:length(f_ex)
    subplot(2,3,f)
    semilogx(amps(freqs==f_ex(f)),n_act(n_ex,freqs==f_ex(f)),'k','LineWidth',1.5)
    hold on
    semilogx(amps(freqs==f_ex(f)),n_opt(n_ex,freqs==f_ex(f)),'Color',[108 150 157]/255,'LineWidth',1.5)
    semilogx(amps(freqs==f_ex(f)),n_test1112(n_ex,freqs==f_ex(f)),'Color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
    box off
    xlim([1 100])
    ylim([0 60])
    xlabel('Amplitude [um]')
    ylabel('Firing rate [Hz]')
    
    if f==1
        legend('data','full model','reduced model')
    end
end

% decay clamped
subplot(2,3,4)
loglog(freqs_unique,thres_act(n_ex,:),'k','LineWidth',1.5)
hold on
loglog(freqs_unique,thres_opt(n_ex,:),'Color',[108 150 157]/255,'LineWidth',1.5)
loglog(freqs_unique,thres_test10(n_ex,:),'Color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
box off
xlabel('Frequency [Hz]')
ylabel('Amplitude [um]')

% no saturation
n_ex = 4;
f_ex = [400 500];
for f=1:length(f_ex)
    subplot(2,3,4+f)
    semilogx(amps(freqs==f_ex(f)),n_act(n_ex,freqs==f_ex(f)),'k','LineWidth',1.5)
    hold on
    semilogx(amps(freqs==f_ex(f)),n_opt(n_ex,freqs==f_ex(f)),'Color',[108 150 157]/255,'LineWidth',1.5)
    semilogx(amps(freqs==f_ex(f)),n_test8(n_ex,freqs==f_ex(f)),'Color',[0.4940 0.1840 0.5560],'LineWidth',1.5)
    box off
    xlim([.1 100])
    ylim([0 80])
    xlabel('Amplitude [um]')
    ylabel('Firing rate [Hz]')
end

screenshot(1,'param_check1','pdf')
