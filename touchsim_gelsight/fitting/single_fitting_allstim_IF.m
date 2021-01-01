% single_fitting_allstim_IF.m

fit_class = 'SA1';
fit_ind = 1;


% SA: 2 4 5  8
% RA: 1 3 11 12 13 14 19 20 22
% PC: 6 7 9  25

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

aff_idx = find(strcmp(aff_class,fit_class(1:2)));
a_ind = aff_idx(fit_ind);

aff_ind_noise = idx(1,a_ind);
aff_ind_noise_lo = idx(2,a_ind);
aff_ind_sine = idx(3,a_ind);
aff_ind_diharm = idx(4,a_ind);

%% set initial parameters

parameters = IF_parameters();

opts = optimset('display', 'iter', 'maxiter', 25, 'maxfunevals', 20000);
%opts = gaoptimset('display', 'iter','InitialPopulation',cell2mat(parameters.ra'));

switch fit_class
    case 'SA1'
        pp_ini = parameters.sa{fit_ind};
        tc = 250;
        %   1:cutoff 2:7:wk 8:sat 9:noise 10:tau 11:12:ih 13:delay
        LB = [2  -10 -10 -1 -1 0 0 0 0 0.01 -Inf -Inf 1];
        UB = [11  10  10  1  0 0 0 0 0 500   0    0   20];
        
    case 'RA'
        pp_ini = parameters.ra{fit_ind};
        tc = 250;
        %   1:cutoff 2:7:wk 8:sat 9:noise 10:tau 11:12:ih 13:delay
        LB = [2  0 0 0   0   -Inf -Inf 0.1  0 0.01 -Inf -Inf 1];
        UB = [20 0 0 Inf Inf  Inf  Inf 1000 0 55    0    0   20];
        
    case 'PC'
        pp_ini = parameters.pc{fit_ind};
        tc = 100;
        %   1:cutoff 2:7:wk 8:sat 9:noise 10:tau 11:12:ih 13:delay
        LB = [2  0 0 -Inf -Inf -Inf -Inf 0.1  0 0.01 -Inf -Inf 1];
        UB = [50 0 0  Inf  Inf  Inf  Inf 1000 0 500   0    0   20];
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
trace = s.propagate(Afferent(fit_class,[0 0]),0);
stim_strain = trace.stat_comp(1:end-2000);
stim_udyn = trace.dyn_comp(1:end-2000);

s = Stimulus(stim2/1000,[0 0],20000,0.5);
trace = s.propagate(Afferent(fit_class,[0 0]),0);
stim_strain = [stim_strain; trace.stat_comp; zeros(2000,1)];
stim_udyn = [stim_udyn; trace.dyn_comp; zeros(2000,1)];
clear s trace

stim(:,1) = reshape(stim_strain,[],1);
stim(:,2) = reshape(stim_udyn,[],1);

% get response
if ~strcmpi(fit_class,'PC')
    [sp,spvec] = get_spikes_mult({NOISE_LO SINE},...
        {aff_ind_noise_lo aff_ind_sine},...
        {idx_noise_lo idx_sine},...
        {repmat([0.5 1.02],length(idx_noise_lo),1) [repmat(0.05,length(idx_sine),1) stim_len+0.02]});
else
    [sp,spvec] = get_spikes_mult({L2_str NOISE_LO SINE},...
        {aff_ind_noise aff_ind_noise_lo aff_ind_sine},...
        {idx_noise idx_noise_lo idx_sine},...
        {repmat([0.5 1.02],length(idx_noise),1) repmat([0.5 1.02],length(idx_noise_lo),1) [repmat(0.05,length(idx_sine),1) stim_len+0.02]});
end

%%

d0 = VRdist_IF(zeros(1,length(pp_ini)),stim,spvec,tc);
fprintf(['Distance with 0 spikes: ' num2str(d0) '\n'])

% optimize parameters
pp = optimize_VRdist_IF(pp_ini,stim,tc,spvec,opts,{},LB,UB);

%%

[pred,predvec] = predict_IF(pp,stim);

fprintf(['Fit complete. ' num2str(length(pred)) ' sim vs ' num2str(length(sp)) ' act spikes.\n'])

% make some figures
figure
plot_fit(stim,spvec,predvec)

figure
hold on
plot(sum(reshape(spvec,75,[])),'r')
plot(sum(reshape(predvec,75,[])),'g')

%% save results

switch fit_class
    case 'SA1'
        parameters.sa{fit_ind} = pp;
    case 'RA'
        parameters.ra{fit_ind} = pp;
    case 'PC'
        parameters.pc{fit_ind} = pp;
end

Simulink.saveVars('../base/internal/IF_parameters.m','parameters')
