function [PC_mean, RA_mean, SA_mean] = plotFiringRates(rates, names, neuron_identities, texture_num, speed, excludeNeurons)
%plotFiringRates plots mean firing rates for each afferent for a given
%texture at a given speed. If exclude neurons is positive, excludes neurons
%that don't fire.

iPC = neuron_identities{1};
iRA = neuron_identities{2};
iSA = neuron_identities{3};
name_str = names{texture_num};

ind = speed/40;
rates_mat = rates{ind};
rates_mat = squeeze(rates_mat(texture_num,:, :));

%alert user to how good the readings are here
percent_non_nan = sum(~isnan(rates_mat(:)))/length(rates_mat);
disp(strcat(num2str(percent_non_nan, "% of neuron responses non-zero")));

stats = {}; % PC RA SA1 are columns. Mean and SEM are rows.
pcs = rates_mat(iPC, :);
pcs_flat = pcs(:);
ras = rates_mat(iRA, :);
ras_flat = ras(:);
sas = rates_mat(iSA, :);
sas_flat = sas(:);

if excludeNeurons
    pcs_flat = pcs_flat(pcs_flat ~= 0);
    ras_flat = ras_flat(ras_flat ~= 0);
    sas_flat = sas_flat(sas_flat ~= 0);
end

stats{1, 1} = nanmean(pcs_flat);
stats{1, 2} = nanmean(ras_flat);
stats{1, 3} = nanmean(sas_flat);

stats{2, 1} = nanstd(pcs_flat)/sqrt(length(~isnan(pcs)));
stats{2, 2} = nanstd(ras_flat)/sqrt(length(~isnan(ras)));
stats{2, 3} = nanstd(sas_flat)/sqrt(length(~isnan(sas)));
x = [1,2,3];

means = [stats{1, 1}, stats{1, 2}, stats{1, 3}];
sem = [stats{2, 1}, stats{2, 2}, stats{2, 3}];
figure;
bar(x, means);
hold on
er = errorbar(x,means,means-sem,means+sem);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title(strcat(name_str, " firing rate, ", num2str(speed), " mm/s"));
xticks(x)
xticklabels({'PCs','RAs','SAs'})
ylabel("Hz")
PC_mean = stats{1,1};
RA_mean = stats{1,2};
SA_mean = stats{1,3};
