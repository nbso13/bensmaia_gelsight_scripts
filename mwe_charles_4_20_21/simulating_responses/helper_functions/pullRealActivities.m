function [activities, spike_trains, space_vec, raw_rates] = pullRealActivities(rates, ...
    spikes, names, good_neurons, neuron_identities, texture_nums, speed, ...
    excludeNeurons, time_samp_period)
%pullRealRates grabs real neural data, and calculated rates for each
%neuron, sorts into afferent class, calculates mean firing rate for each
%texture for each afferent, and returns a mean_rates cell array, where rows are
%textures and columns are arrays of afferent classes (PC, RA, SA) and entries are
%average firing rates for each neuron. Next three columns are arrays of standard error on the mean of
%same afferent classes.

% for spike_trains, returns a cell array where rows are textures
% and columns are afferent classes. Entries are cell arrays where each cell
% is a neuron. They themselves contain cell arrays where each array is a
% a successful trial where entries are impulses rather than spike times, x axis
% being in space. This is for later correlation with stimulus PSD.

%output raw rates is a cell array, texture_num by 3. Each entry is a vector
%of rates for that aff and that texture.

if texture_nums == 1
    error("for now, this doesn't work with only one texture.")
end

%determine "good neuron" indices
good_logit = zeros(39,1);
for i = 1:length(good_logit)
    if ismember(i, good_neurons)
        good_logit(i) = 1;
    end
end
good_logit = logical(good_logit);

%indices for each class
iPC = neuron_identities{1};
iRA = neuron_identities{2};
iSA = neuron_identities{3};

iPC = iPC(good_logit); % choosing good neurons
iRA = iRA(good_logit);
iSA = iSA(good_logit);

ind = speed/40;
spikes_mat = spikes{ind};
spikes_mat = squeeze(spikes_mat(texture_nums,:, :));
spikes_mat = spikes_mat(:, good_logit, :);


pcs = spikes_mat(:, iPC, :);
ras = spikes_mat(:, iRA, :);
sas = spikes_mat(:, iSA, :);
afferents = {pcs, ras, sas};

samp_period = time_samp_period; %seconds LATER MULT BY 80mm/s
dec_place = -1*log10(samp_period);
time_axis = 0+samp_period:samp_period:1;

spike_trains = {};

for i = 1:length(texture_nums) %for every texture
    for j = 1:3 %for each different afferent type
        aff = afferents{j};
        aff_array = zeros(1, length(time_axis));
        for k = 1:size(aff, 2) %for every neuron of that afferent
            for m = 1:4 %for each possible trial
                run = aff{i,k,m}; %grab the spike times for that texture, aff type, and trial
                if ~isempty(run) %if this is not empty
                    run = round(run, dec_place); %round spike times as appropriate
                    run = run.*(10^dec_place)+1; %turning spike times into indices
                    run_vector = zeros(1, length(time_axis));
                    run_vector(run) = 1; %add one at each appropriate index
                    run_vector = run_vector(1:length(time_axis));
                    aff_array = vertcat(aff_array, run_vector);
                end
                
            end
            
        end
        spike_trains{i, j} = aff_array(2:end, :); %to not include initial row of 0's
    end
end

space_vec = time_axis.*speed; %now x vector is in mm.

%% REPEATING FOR RATES

rates_mat = rates{ind};
rates_mat = squeeze(rates_mat(texture_nums,:, :));
rates_mat = rates_mat(:, good_logit,:);
raw_rates = {};

% nanmean trials together for each neuron for each texture.
rates_mat_meaned = nanmean(rates_mat, 3);



% rate_mat is now two dims - texture num, neuron num. Entries are mean rates
% over 4 trials for that neuron for that texture. Here we're
% collapsing the last two dimensions to give a cell array, each cell an
% afferent type, each cell texture_num by number of afferents.
full_rates = {reshape(rates_mat_meaned(:, iPC, :), [length(texture_nums), sum(iPC)]),...
    reshape(rates_mat_meaned(:, iRA, :), [length(texture_nums), sum(iRA)]),...
    reshape(rates_mat_meaned(:, iSA, :), [length(texture_nums), sum(iSA)])};

for i = 1:length(full_rates) %for each afferent type
    aff_rates =  full_rates{i};
    if excludeNeurons
        aff_rates(aff_rates == 0) = nan;
    end
    mean_rates(:, i) = mean(aff_rates, 2, 'omitnan');
    mean_rates(:, i+3) = std(aff_rates, 0, 2, 'omitnan')/sqrt(length(~isnan(aff_rates)));
    for j = 1:length(texture_nums) 
        %setting up raw rates, cell array texture_num by 3, each entry a vector of non-nan rates for that aff and that texture
         aff_texture_rates = aff_rates(j, :);
%          aff_texture_rates = aff_texture_rates(logical(~isnan(aff_texture_rates)));
         %grab non-nan
         raw_rates{j, i} = aff_texture_rates;
    end
end

activities = struct;
activities.names = names;
activities.real = mean_rates;
end

