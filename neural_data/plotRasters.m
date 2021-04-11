function [] = plotRasters(spikes, time_length, good_neurons, names, neuron_identities, texture_num, speed)
%plotRasters takes spiketrains for a given texture at a given speed and
%plots the raster. SAs are green, RAs are blue, PCs are orange.
good_logit = zeros(39,1);
for i = 1:length(good_logit)
    if ismember(i, good_neurons)
        good_logit(i) = 1;
    end
end
good_logit = logical(good_logit);

iPC = neuron_identities{1};
iRA = neuron_identities{2};
iSA = neuron_identities{3};

iPC = iPC(good_logit);
iRA = iRA(good_logit);
iSA = iSA(good_logit);
name_str = names{texture_num};

ind = speed/40;
spikes_mat = spikes{ind};
spikes_mat = squeeze(spikes_mat(texture_num,:, :));
spikes_mat = spikes_mat(good_logit, :);

pcs = spikes_mat(iPC, :);
ras = spikes_mat(iRA, :);
sas = spikes_mat(iSA, :);
neurons = {sas, ras, pcs};
colors = [[0.2 0.85 0.2]; [0.2 0.2 1]; [1 0.62 0]];
hold on
plotter = 1;
for i = 1:length(neurons)
    type = neurons{i};
    for j = 1:size(type, 1)
        for k = 1:size(type, 2)
            spike_train = type{j,k};
            if size(spike_train, 1) == 0
                continue
            else
                for m = 1:length(spike_train)
                    if spike_train(m) > time_length
                        break
                    end
                    line([spike_train(m) spike_train(m)], [plotter-1 plotter], 'Color', colors(i, :));
                    %Create a tick mark at x = t1(i) with a height of 1
                    
                    
                end
                plotter = plotter + 1;
            end
        end
    end
    plotter = plotter + 4;
            
xlabel('Time (sec)') %Label x-axis
ylabel('Neurons')

title(strcat(name_str, " Raster, ", num2str(speed), " mm/s"));
end

