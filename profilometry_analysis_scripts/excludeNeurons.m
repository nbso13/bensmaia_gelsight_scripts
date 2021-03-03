function [new_res_coll, new_affpop] = excludeNeurons(res, top_neurons, aff_pop)
%excludeNeurons takes only the responses from top n neurons in each
%afferent
% top_neurons_number = top_neurons_max - top_neurons_min;
len_rates = length(res.rate);
indices = 1:len_rates;
rate_array = {res.rate(aff_pop.iPC), res.rate(aff_pop.iRA),res.rate(aff_pop.iSA1)}; % separate list of rates
master_indices_array = {indices(aff_pop.iPC), indices(aff_pop.iRA), indices(aff_pop.iSA1)}; % separate lists of indices
master_top_indices_list = [];
for i = 1:length(rate_array)
    [~, sorted_indices] = sort(rate_array{i}, 'descend');
    % 10 median neurons
    top_indices = sorted_indices(1:top_neurons); %top indices IN EACH AFFERENT CATEGORY
    number_total_afferents = length(sorted_indices);
    aff_indices_mask = zeros(1, number_total_afferents); % make mask to find the top master indices
    aff_indices_mask(top_indices) = 1; %set mask to one in found top indices
    aff_indices_mask = logical(aff_indices_mask);
    master_aff_ind_arr = master_indices_array{i}; % get master indices list for this afferent
    top_master_indices = master_aff_ind_arr(aff_indices_mask); %take according to local mask
    master_top_indices_list = [master_top_indices_list, top_master_indices]; %append indices onto master list of top indices
end

%make overall mask
overall_top_mask = zeros(1,len_rates);
overall_top_mask(master_top_indices_list) = 1;
overall_top_mask = logical(overall_top_mask);

%restricting responses
new_responses = res.responses(overall_top_mask);



%restricting affpop
new_affpop = res.affpop;
new_affpop.afferents = new_affpop.afferents(overall_top_mask);
% new_affpop.location = new_affpop.location(good_rate);
% new_affpop.depth = new_affpop.depth(good_rate);
% new_affpop.num = length(good_rate);
% new_affpop.class = new_affpop.class(good_rate);
% new_affpop.iSA1 = new_affpop.iSA1(good_rate);
% new_affpop.iRA = new_affpop.iRA(good_rate);
% new_affpop.iPC = new_affpop.iPC(good_rate);




new_res_coll = ResponseCollection(new_affpop, new_responses, res.stimulus);
end

