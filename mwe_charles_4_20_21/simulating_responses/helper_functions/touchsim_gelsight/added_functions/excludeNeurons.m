function [new_res_coll] = excludeNeurons(res,rate_lower_lim)
%excludeNeurons takes only the responses with rates equal to or higher than
%rate lower lim. discards

%get indices of neurons that fire at or above rate
good_rate = res.rate>=rate_lower_lim;

%restricting responses
new_responses = res.responses(good_rate);

%restricting rate
new_rate = res.rate(good_rate);



%restricting affpop
new_affpop = res.affpop;
new_affpop.afferents = new_affpop.afferents(good_rate);
% new_affpop.location = new_affpop.location(good_rate);
% new_affpop.depth = new_affpop.depth(good_rate);
% new_affpop.num = length(good_rate);
% new_affpop.class = new_affpop.class(good_rate);
% new_affpop.iSA1 = new_affpop.iSA1(good_rate);
% new_affpop.iRA = new_affpop.iRA(good_rate);
% new_affpop.iPC = new_affpop.iPC(good_rate);




new_res_coll = ResponseCollection(new_affpop, new_responses, res.stimulus);
end

