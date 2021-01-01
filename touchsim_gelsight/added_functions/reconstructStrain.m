function [strain_at_time] = reconstructStrain(r, time_point)
%ReconstructStrain This func visualizes the local strain for each receptor
%based on encoded location on the finger. Time point in seconds.

dur = r.responses(1,1).duration;
samp_freq = r.responses(1,1).propagated_struct.sampling_frequency;
num_neurons = length(r.responses);

%first we need the total trace for each neuron
total_trace_mat = zeros(num_neurons, samp_freq*dur);
for i = 1:num_neurons
    total_trace_mat(i, :) = r.responses(1,i).propagated_struct.stat_comp(:);
end

strain_at_time = total_trace_mat(:,time_point*samp_freq);

% 
% %dimesions of strainmat are the max in both dimensions of the locations of
% %the afferents.
% locations_scaled = floor(r.affpop.location.*100);
% locations_demin = locations_scaled - min(locations_scaled, [], 'all')+1;
% y_max = max(locations_demin(:,1), [], 'all');
% x_max = max(locations_demin(:,2), [], 'all');
% strainmat = zeros(y_max, x_max);
% 
% %assign strain value at time point
% for i = 1:num_neurons
%     x = locations_demin(i, 1);
%     y = locations_demin(i, 2);
%     strainmat(x,y) = strain_at_time(i);
% end


end

