%% Analyzing Neural Responses Periphery to Textures

clear 
close all

load('RawPAFData.mat')
load('TextureNames')

% About the data -
% good neurons indices [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 22 25 28 33 34]
% Recorded from 39 neurons. 55 textures were spun across monkey fingers.
% Four reps at three speeds (40, 80, 120 mm/s I believe)
% good_neurons = [2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 18 22 25 28 33 34];
good_neurons = 1:39;
% text_nums = [6, 9, 13, 21, 25, 37, 31, 50:55];

text_nums = [7, 31, 4, 42, 45, 50, 49, 48];
sub_plot_y_len = 3;
sub_plot_x_len = 3;


activities = struct;

activities.names = ["wool_blend","velvet","hucktowel","sueded_cuddle","blizzard_fleece","1mm_grating","3mm_grating", "5mm_grating"];
activities.real = zeros(length(activities.names), 6);
speeds = [80];
for j = 1:length(speeds)
    speed  = speeds(j);
    neuron_identities = {iPC, iRA, iSA};
    figure('Position', [100, 100, 500, 500]);
    for i = 1:length(text_nums)
        subplot(sub_plot_y_len,sub_plot_x_len,i)
        [means, sems] = plotFiringRates(rates, good_neurons, htxt_name, neuron_identities, text_nums(i), speed, 1);
        activities.real(i,:) = [means, sems];
        ylim([0 inf])
        
    end
    
    figure('Position', [400, 400, 800, 500]);
    for i = 1:length(text_nums)
        subplot(sub_plot_y_len,sub_plot_x_len,i)
        plotRasters(spikes, good_neurons, htxt_name, neuron_identities, text_nums(i), speed);
    end
end


ind = speed/40;
rates_mat = rates{ind};
rates_mat = squeeze(rates_mat(texture_num,:, :));
stats = {}; % PC RA SA1 are columns. Mean and SEM are rows.
pcs = rates_mat(iPC, :);
ras = rates_mat(iRA, :);
sas = rates_mat(iSA, :);

stats{1, 1} = nanmean(pcs(:));
stats{1, 2} = nanmean(ras(:));
stats{1, 3} = nanmean(sas(:));

stats{2, 1} = nanstd(pcs(:))/sqrt(length(~isnan(pcs)));
stats{2, 2} = nanstd(ras(:))/sqrt(length(~isnan(ras)));
stats{2, 3} = nanstd(sas(:))/sqrt(length(~isnan(sas)));

figure
x = [1,2,3];

means = [stats{1, 1}, stats{1, 2}, stats{1, 3}];
sem = [stats{2, 1}, stats{2, 2}, stats{2, 3}];
figure;
bar(x, means);
hold on
er = errorbar(x,means,means-sem,means+sem);
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
title(strcat(htxt_name{texture_num}, " firing rate"));
xticks(x)
xticklabels({'PCs','RAs','SAs'})
ylabel("Hz")




for i = 1:2
    disp(htxt_name(texture_num))
    temp = spikes{1}(texture_num,num,i);
    spike_train_last_texture{1,i} = temp;
end

% grab rates for 2 mm dot patterns for good neurons at 80mm/s
data = [];
twommdots_eightymms = {};
data_spikes = squeeze(spikes{2}(texture_num,:,:)); %80m/s, only for dots
data_rate = squeeze(rates{2}(texture_num,:,:));
good_SAs = {};
good_SA_count = 1;
for k = 1:length(indices)
    i = indices(k); %i is the index of a good neuron
    if iSA(i) %if an SA
        neuron_count = 1;
        for j = 1:4 %for all runs
            temp_var = data_spikes(i, j);
            temp_var = temp_var{1};
            if ~isempty(temp_var) % if there are spikes here
                temp{neuron_count} = temp_var;  %put spikes in temp
                temp_rates{neuron_count} = data_rate(i, j); %put rates in temp
                neuron_count= neuron_count+1;
            end
        end
        if neuron_count > 1 %if there was at least one good trial
            for j = 1:length(temp)
                good_SAs{good_SA_count, j} = temp{j};
                SA_rates(good_SA_count, j) = temp_rates{j};
            end
            good_SA_count = good_SA_count+1;
        end
    end
    
end
% now good_SAs and SA_rates are organized like this  - each row vector is a
% neuron and each column is a trial. If a neuron has no activity at all it
% is left out. No trials with no activity are included.

figure
plotter_count = 1;
hold on %Allows multiple plots on the same graph
for i = 1:size(good_SAs, 1) % for each neuron
    for k = 1:size(good_SAs, 2)
        temp = good_SAs{i, k};
        for j = 1:length(temp) %Going through each spike time
            line([temp(j) temp(j)], [plotter_count-1 plotter_count]); %Create a tick mark at x = t1(i) with a height of 1
        end
        plotter_count = plotter_count+1;
    end
    plotter_count=plotter_count+1;
    
end

ylim([0 31]) %Reformat y-axis for legibility
xlabel('Time (sec)') %Label x-axis
ylabel('Trial #')

title(strcat("Good SAs for ", htxt_name{texture_num}, " Scanned at 80mm/s (Weber et al 2013)"))


for i = 1:length(SA_rates)
    mean_fr(i) = mean(SA_rates(i,:));
    sd_fr(i) = std(SA_rates(i,:));
end

figure

x = 1:length(mean_fr);
errhigh = sd_fr;
errlow = sd_fr;

bar(x, mean_fr)

hold on

er = errorbar(x, mean_fr, errlow, errhigh);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
title(strcat('Mean FR in Hz Across Trials for Each Neuron (Tx ', num2str(texture_num),  ' sp=80)'))
xlabel('Neuron')
ylabel('FR (Hz)')

figure
m = mean(mean_fr);
s = std(mean_fr);
bar(1, m)
hold on
er = errorbar(1, m, s, s);
er.Color = [0 0 0];                            
er.LineStyle = 'none'; 
title(strcat('Mean FR in Hz (Tx ', num2str(texture_num),' sp=80)'))
xlabel('population average')
ylabel('FR (Hz)')


WRITE ANALYSIS CODE HERE