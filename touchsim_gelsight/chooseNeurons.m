function [resp_new, aff_pop_new] = chooseNeurons(resp, neuron_selection_modes, ...
        texture_rates, top_neuron_number, aff_pop, loc)
%chooseNeurons selects neurons in the population based on their location or responses.

min_x = loc(1);
max_x = loc(2);
min_y = loc(3);
max_y = loc(4);

for i = 1:3 %for each afferent
    if strcmp(neuron_selection_modes(i), "top")
        [resp_new, aff_pop_new] = excludeNeurons(resp, top_neuron_number, aff_pop, i);
        
        
    elseif strcmp(neuron_selection_modes(i), "area")
        
        
    elseif strcmp(neuron_selection_modes(i), "best")
        
        
    else
        error("neuron selection mode not recognized")
    end

end
end

