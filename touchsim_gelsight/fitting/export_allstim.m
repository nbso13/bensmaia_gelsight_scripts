% export_allstim_IF.m

load('pma_noise_chub.mat');
NOISE_LO = L2_str;
load('pma_chubNoise_desiredpositions.mat')
load('pma_noise.mat')
load('SINE.mat')
load('DI.mat')

[idx,aff_class] = find_common_neurons(L2_str,NOISE_LO,SINE,DI);

sub_ind = ~isnan(sum(idx));

idx = idx(:,sub_ind);
aff_class = aff_class(sub_ind);

aff_ind_noise = idx(1,:);
aff_ind_noise_lo = idx(2,:);
aff_ind_sine = idx(3,:);
aff_ind_diharm = idx(4,:);

%%

indSA = strmatch('SA',aff_class);
indRA = strmatch('RA',aff_class);
indPC = strmatch('PC',aff_class);

% extract stimulus
stim_noise = cell2mat(L2_str.stim_pos{1}')';
idx_noise = 1:40;

stim_noise_lo = cell2mat(desiredpositions');
stim_noise_lo = stim_noise_lo(:,10001:30000);
idx_noise_lo = 1:20;

stim_len = SINE.state_v(:,8);
idx_sine = 1:length(stim_len);
[freqs,ind_sorted] = sort(SINE.state_v(idx_sine,5),'ascend');
idx_sine = idx_sine(ind_sorted);
stim_len = stim_len(ind_sorted);

load desiredsinepos.mat
for ii=1:length(desiredsinepos)
    [b,a]=butter(1,min(freqs(ii)*8/10000,.99));
    stim_sine_all{ii}=filter(b,a,[zeros(1,2000) desiredsinepos{ii} zeros(1,2000)]);
end

idx_sineRA = idx_sine(1:120);
stim_lenRA = stim_len(1:120);
stim_sine_allRA = stim_sine_all(1:120);

stim1 = get_stimulus(4000,stim_noise,stim_noise_lo);
s = Stimulus(stim1/1000,[0 0],20000,0.5);
stim1RA = get_stimulus(4000,stim_noise_lo);
sRA = Stimulus(stim1RA/1000,[0 0],20000,0.5);
trace = sRA.propagate(Afferent('SA1',[0 0]),0);
stimSA1 = trace.stat_comp(1:end-2000);
stimSA2 = trace.dyn_comp(1:end-2000);
trace = sRA.propagate(Afferent('RA',[0 0]),0);
stimRA1 = trace.stat_comp(1:end-2000);
stimRA2 = trace.dyn_comp(1:end-2000);
trace = s.propagate(Afferent('PC',[0 0]),0);
stimPC1 = trace.stat_comp(1:end-2000);
stimPC2 = trace.dyn_comp(1:end-2000);

stim2 = get_stimulus(0,stim_sine_all{:});
s = Stimulus(stim2/1000,[0 0],20000,0.5);
stim2RA = get_stimulus(0,stim_sine_allRA{:});
sRA = Stimulus(stim2RA/1000,[0 0],20000,0.5);
trace = sRA.propagate(Afferent('SA1',[0 0]),0);
stimSA1 = [stimSA1; trace.stat_comp; zeros(2000,1)];
stimSA2 = [stimSA2; trace.dyn_comp; zeros(2000,1)];
trace = sRA.propagate(Afferent('RA',[0 0]),0);
stimRA1 = [stimRA1; trace.stat_comp; zeros(2000,1)];
stimRA2 = [stimRA2; trace.dyn_comp; zeros(2000,1)];
trace = s.propagate(Afferent('PC',[0 0]),0);
stimPC1 = [stimPC1; trace.stat_comp; zeros(2000,1)];
stimPC2 = [stimPC2; trace.dyn_comp; zeros(2000,1)];
clear s trace

stimi = decimate(reshape(stimSA1,[],1),4);
dstimi = decimate(reshape(stimSA2,[],1),4);
ddstimi = [diff(dstimi);0];
inputSA = rectify([stimi -stimi dstimi -dstimi ddstimi -ddstimi]);

stimi = decimate(reshape(stimRA1,[],1),4);
dstimi = decimate(reshape(stimRA2,[],1),4);
ddstimi = [diff(dstimi);0];
inputRA = rectify([stimi -stimi dstimi -dstimi ddstimi -ddstimi]);

stimi = decimate(reshape(stimPC1,[],1),4);
dstimi = decimate(reshape(stimPC2,[],1),4);
ddstimi = [diff(dstimi);0];
inputPC = rectify([stimi -stimi dstimi -dstimi ddstimi -ddstimi]);

for a=1:length(aff_ind_noise)
    
    % get response
    if ~strcmpi(aff_class{a},'PC')
        [sp{a},spvec{a}] = get_spikes_mult({NOISE_LO SINE},...
            {aff_ind_noise_lo(a) aff_ind_sine(a)},...
            {idx_noise_lo idx_sineRA},...
            {repmat([0.5 1.02],length(idx_noise_lo),1) [repmat(0.05,length(idx_sineRA),1) stim_lenRA+0.02]});
    else
        [sp{a},spvec{a}] = get_spikes_mult({L2_str NOISE_LO SINE},...
            {aff_ind_noise(a) aff_ind_noise_lo(a) aff_ind_sine(a)},...
            {idx_noise idx_noise_lo idx_sine},...
            {repmat([0.5 1.02],length(idx_noise),1) repmat([0.5 1.02],length(idx_noise_lo),1) [repmat(0.05,length(idx_sine),1) stim_len+0.02]});
    end
end

%% save results

save training_data.mat aff_class spvec sp inputSA inputRA inputPC
