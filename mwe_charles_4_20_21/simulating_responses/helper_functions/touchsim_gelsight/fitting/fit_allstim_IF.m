% fit_allstim_IF.m

fit_class = '';

NOISE_LO = load('e:\Dropbox\chi\Matlab scripts\peripheral\mat\pma_noise_chub.mat');
NOISE_LO = NOISE_LO.L2_str;
load('e:\Dropbox\chi\Matlab scripts\peripheral\mat\pma_chubNoise_desiredpositions.mat')
load('e:\Dropbox\chi\Matlab scripts\peripheral\mat\pma_noise.mat')
load('e:\Dropbox\chi\Matlab scripts\chicago2\mat\SINE.mat')
load('e:\Dropbox\chi\Matlab scripts\chicago2\mat\DI.mat')

[idx,aff_class] = find_common_neurons(L2_str,NOISE_LO,SINE,DI);

if isempty(fit_class)
    sub_ind = ~isnan(sum(idx));
else
    sub_ind = ~isnan(sum(idx)) & strcmpi(fit_class,aff_class)';
end
idx = idx(:,sub_ind);
aff_class = aff_class(sub_ind);

aff_ind_noise = idx(1,:);
aff_ind_noise_lo = idx(2,:);
aff_ind_sine = idx(3,:);
aff_ind_diharm = idx(4,:);

%% set initial parameters

IF_parameters;

opts = optimset('display', 'iter', 'maxiter', 150, 'maxfunevals', 20000);
%opts = gaoptimset('display', 'iter','PopulationSize',21,'InitialPopulation',[cell2mat(parameters.sa'); cell2mat(parameters.ra'); cell2mat(parameters.pc')]);
%opts = gaoptimset('display', 'iter','InitialPopulation',[cell2mat(parameters.sa'); cell2mat(parameters.ra'); cell2mat(parameters.pc')]);

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
trace = sRA.propagate(Afferent('SA1',[0 0]));
stimSA1 = trace.stat_comp(1:end-2000);
stimSA2 = trace.dyn_comp(1:end-2000);
trace = sRA.propagate(Afferent('RA',[0 0]));
stimRA1 = trace.stat_comp(1:end-2000);
stimRA2 = trace.dyn_comp(1:end-2000);
trace = s.propagate(Afferent('PC',[0 0]));
stimPC1 = trace.stat_comp(1:end-2000);
stimPC2 = trace.dyn_comp(1:end-2000);

stim2 = get_stimulus(0,stim_sine_all{:});
s = Stimulus(stim2/1000,[0 0],20000,0.5);
stim2RA = get_stimulus(0,stim_sine_allRA{:});
sRA = Stimulus(stim2RA/1000,[0 0],20000,0.5);
trace = sRA.propagate(Afferent('SA1',[0 0]));
stimSA1 = [stimSA1; trace.stat_comp; zeros(2000,1)];
stimSA2 = [stimSA2; trace.dyn_comp; zeros(2000,1)];
trace = sRA.propagate(Afferent('RA',[0 0]));
stimRA1 = [stimRA1; trace.stat_comp; zeros(2000,1)];
stimRA2 = [stimRA2; trace.dyn_comp; zeros(2000,1)];
trace = s.propagate(Afferent('PC',[0 0]));
stimPC1 = [stimPC1; trace.stat_comp; zeros(2000,1)];
stimPC2 = [stimPC2; trace.dyn_comp; zeros(2000,1)];
clear s trace


stimSA(:,1) = reshape(stimSA1,[],1);
stimSA(:,2) = reshape(stimSA2,[],1);
stimRA(:,1) = reshape(stimRA1,[],1);
stimRA(:,2) = reshape(stimRA2,[],1);
stimPC(:,1) = reshape(stimPC1,[],1);
stimPC(:,2) = reshape(stimPC2,[],1);

tic
pp = cell(length(aff_ind_noise),1);
ppa = cell(length(aff_ind_noise),1);
pred = cell(length(aff_ind_noise),1);
pred_vec = cell(length(aff_ind_noise),1);

cluster = parpool(7);
parfor a=1:length(aff_ind_noise)
%for a=1:length(aff_ind_noise)
    
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
    
    
    switch aff_class{a}
        case 'SA'
            idx = find(indSA==a);
            ppa{a} = parameters.sa{idx};
        case 'RA'
            idx = find(indRA==a);
            ppa{a} = parameters.ra{idx};
        case 'PC'
            idx = find(indPC==a);
            ppa{a} = parameters.pc{idx};
    end
    
    % optimize parameters
    fprintf(['Fit ' num2str(a) ' (' aff_class{a} ') starting.\n'])
    switch aff_class{a}
        case 'SA'
            LB = [0.01 -1 -1 -1 -1 0 0 -Inf -Inf 0  1  2 ];
            UB = [500   1  1  1  1 0 0  0    0   0 20  11];
            %pp{a} = optimize_VRdist_IF(ppa{a},stimSA,1000,spvec{a},opts,parameters.sa,LB,UB);
            pp{a} = optimize_VRdist_IF(ppa{a},stimSA,250,spvec{a},opts,{},LB,UB);
            [pred{a},predvec{a}] = predict_IF(pp{a},stimSA);
        case 'RA'
            LB = [0.01 0 0 -Inf -Inf -Inf -Inf -Inf -Inf 0.1  1  2 ];
            UB = [500  0 0  Inf  Inf  Inf  Inf  0    0   1000 20 20];
            %pp{a} = optimize_VRdist_IF(ppa{a},stimRA,1000,spvec{a},opts,parameters.ra,LB,UB);
            pp{a} = optimize_VRdist_IF(ppa{a},stimRA,250,spvec{a},opts,{},LB,UB);
            [pred{a},predvec{a}] = predict_IF(pp{a},stimRA);
        case 'PC'
            LB = [0.01 0 0 -Inf -Inf -Inf -Inf -Inf -Inf 0.1  1  2 ];
            UB = [500  0 0  Inf  Inf  Inf  Inf  0    0   1000 20 50];
            %pp{a} = optimize_VRdist_IF(ppa{a},stimPC,100,spvec{a},opts,parameters.pc,LB,UB);
            pp{a} = optimize_VRdist_IF(ppa{a},stimPC,100,spvec{a},opts,{},LB,UB);
            [pred{a},predvec{a}] = predict_IF(pp{a},stimPC);
    end
    
    fprintf(['Fit ' num2str(a) ' (' aff_class{a} ') complete. ' num2str(length(pred{a})) ' sim vs ' num2str(length(sp{a})) ' act spikes.\n'])
end
delete(cluster)
toc

%% collate and save models

if isempty(fit_class) || strcmpi(fit_class,'SA')
    for n=1:length(indSA)
        parameters.sa{n} = pp{indSA(n)};
    end
end

if isempty(fit_class) || strcmpi(fit_class,'RA')
    for n=1:length(indRA)
        parameters.ra{n} = pp{indRA(n)};
    end
end

if isempty(fit_class) || strcmpi(fit_class,'PC')
    for n=1:length(indPC)
        parameters.pc{n} = pp{indPC(n)};
    end
end

parameters.description = {...
    '1: Time constant tau',...
    '2: Positive strain coefficient',...
    '3: Negative strain coefficient',...
    '4: Positive udyn coefficient',...
    '5: Negative udyn coefficient',...
    '6: Positive udyn derivative coefficient',...
    '7: Negative udyn derivative coefficient',...
    '8: 1st coefficient for post-spike current',...
    '9: 2nd coefficient for post-spike current',...
    '10: Saturation parameter',...
    '11: Delay in ms',...
    '12: Low-pass filter cutoff',...
    '13: Noise term'};

Simulink.saveVars('../base/internal/IF_parameters.m','parameters')
