% fit_noise_term.m

fit_class = 'RA';
fit_ind = 9;


% SA: 2 4 5  8
% RA: 1 3 10 11 12 13 14 16 17 18 19 20 21
% PC: 6 7 9  23

NOISE_LO = load('e:\Dropbox\chi\Matlab scripts\peripheral\mat\pma_noise_chub.mat');
NOISE_LO = NOISE_LO.L2_str;
load('e:\Dropbox\chi\Matlab scripts\peripheral\mat\pma_chubNoise_desiredpositions.mat')
load('e:\Dropbox\chi\Matlab scripts\peripheral\mat\pma_noise.mat')
load('e:\Dropbox\chi\Matlab scripts\chicago2\mat\SINE.mat')
load('e:\Dropbox\chi\Matlab scripts\chicago2\mat\DI.mat')

[idx,aff_class] = find_common_neurons(L2_str,NOISE_LO,SINE,DI);

sub_ind = ~isnan(sum(idx));
idx = idx(:,sub_ind);
aff_class = aff_class(sub_ind);

aff_idx = find(strcmp(aff_class,fit_class(1:2)));
a_ind = aff_idx(fit_ind);

aff_ind_noise = idx(1,a_ind);
aff_ind_noise_lo = idx(2,a_ind);
aff_ind_sine = idx(3,a_ind);
aff_ind_diharm = idx(4,a_ind);

%% set initial parameters

IF_parameters;

switch fit_class
    case 'SA1'
        pp = parameters.sa{fit_ind};
        tc = 250;        
    case 'RA'
        pp = parameters.ra{fit_ind};
        tc = 250;
    case 'PC'
        pp = parameters.pc{fit_ind};
        tc = 100;
end

%%

% extract stimulus
stim_noise = cell2mat(L2_str.stim_pos{aff_ind_noise}')';
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

if ~strcmpi(fit_class,'PC')
    idx_sine = idx_sine(1:120);
    stim_len = stim_len(1:120);
    stim_sine_all = stim_sine_all(1:120);
end

if strcmpi(fit_class,'PC')
    stim1 = get_stimulus(4000,stim_noise,stim_noise_lo);    
else
    stim1 = get_stimulus(4000,stim_noise_lo);
end

stim2 = get_stimulus(0,stim_sine_all{:});
stim_strain = [];
stim_udyn = [];

s = Stimulus(stim1/1000,[0 0],20000,0.5);
trace = s.propagate(Afferent(fit_class,[0 0]));
stim_strain = trace.stat_comp(1:end-2000);
stim_udyn = trace.dyn_comp(1:end-2000);

s = Stimulus(stim2/1000,[0 0],20000,0.5);
trace = s.propagate(Afferent(fit_class,[0 0]));
stim_strain = [stim_strain; trace.stat_comp; zeros(2000,1)];
stim_udyn = [stim_udyn; trace.dyn_comp; zeros(2000,1)];
clear s trace

stim(:,1) = reshape(stim_strain,[],1);
stim(:,2) = reshape(stim_udyn,[],1);

% get response
for rep=1:5
    if ~strcmpi(fit_class,'PC')
        [~,spvec{rep}] = get_spikes_mult({NOISE_LO SINE},...
            {aff_ind_noise_lo aff_ind_sine},...
            {idx_noise_lo idx_sine},...
            {repmat([0.5 1.02],length(idx_noise_lo),1) [repmat(0.05,length(idx_sine),1) stim_len+0.02]},rep);
    else
        [~,spvec{rep}] = get_spikes_mult({L2_str NOISE_LO SINE},...
            {aff_ind_noise aff_ind_noise_lo aff_ind_sine},...
            {idx_noise idx_noise_lo idx_sine},...
            {repmat([0.5 1.02],length(idx_noise),1) repmat([0.5 1.02],length(idx_noise_lo),1) [repmat(0.05,length(idx_sine),1) stim_len+0.02]},rep);
    end
end


%% calculate pairwise distances

c = 1;
for r1=1:5
    for r2=(r1+1):5
        dist(c) = VRdist_pairwise(spvec{r1},spvec{r2},tc);
        c = c + 1;
    end
end

%%

[~,pred_orig] = predict_IF(pp,stim);

pp(13) = 0.05;
for rep=1:5
    [~,predvec{rep}] = predict_IF(pp,stim,true);
end

c = 1;
for r1=1:5
    for r2=(r1+1):5
        distpred(c) = VRdist_pairwise(predvec{r1},predvec{r2},tc);
        c = c + 1;
    end
end

bounds = [0 pp(13); 0 mean(distpred)];

for i=1:25;
    fprintf('.')
    noise = mean(bounds(1,:));
    pp(13) = noise;
    
    for rep=1:5
        [~,predvec{rep}] = predict_IF(pp,stim,true);
    end
    c = 1;
    for r1=1:5
        for r2=(r1+1):5
            distpred(c) = VRdist_pairwise(predvec{r1},predvec{r2},tc);
            c = c + 1;
        end
    end
    
    if mean(distpred)>mean(dist)
        bounds(1,2) = noise;
        bounds(2,2) = mean(distpred);
    else
        bounds(1,1) = noise;
        bounds(2,1) = mean(distpred);
    end
end
fprintf('\n')
fprintf(['Final noise term: ' num2str(noise) ', final distance: ', num2str(mean(distpred)) ', target distance: ' num2str(mean(dist)) '\n'])
fprintf(['New model spike count: ' num2str(round(mean(cellfun(@sum,predvec)))) ', Noiseless model spike count: ' num2str(sum(pred_orig))  ', Actual spike count: ' num2str(round(mean(cellfun(@sum,spvec)))) '\n'])

%% generate rasters

figure
hold on
for rep=1:5
    stimes = find(spvec{rep})';
    for ss=1:length(stimes)
        plot([stimes(ss) stimes(ss)],[rep rep+0.75],'Color','k')
    end
end
stimes = find(pred_orig)';
    for ss=1:length(stimes)
        plot([stimes(ss) stimes(ss)],[-0.375 0.375],'Color','b')
    end
for rep=1:5
    stimes = find(predvec{rep})';
    for ss=1:length(stimes)
        plot([stimes(ss) stimes(ss)],[-rep -rep-0.75],'Color','r')
    end
end

%%

switch fit_class
    case 'SA1'
        parameters.sa{fit_ind} = pp;
    case 'RA'
        parameters.ra{fit_ind} = pp;
    case 'PC'
        parameters.pc{fit_ind} = pp;
end

Simulink.saveVars('../base/internal/IF_parameters.m','parameters')
