function [params_new,fval] = optimize_VRdist_IF(pp,stim,tc,spvec,opts,pp_alt,LB,UB)
% [ppnew,prs] = optimize_VRdist(pp,stim_expanded,opts,LB,UB)

fval_old = VRdist_IF(pp,stim,spvec,tc);

pp_all = {pp pp_alt{:}};

for p=1:length(pp_all)

    fprintf(['Run ' num2str(p) ' of ' num2str(length(pp_all)) '.\n'])
    
    % Search for ML estimate
    if nargin==8
        [prs{p},fval(p)] = patternsearch(@(x) VRdist_IF(x,stim,spvec,tc),...
            pp_all{p},[],[],[],[],LB,UB,[],opts);
        %[prs{p},fval(p)] = ga(@(x) VRdist_IF(x,stim,spvec,tc),...
        %    length(pp_all{p}),[],[],[],[],LB,UB,[],opts);
    else
        [prs{p},fval(p)] = patternsearch(@(x) VRdist_IF(x,stim,spvec,tc),...
            pp_all{p},[],[],[],[],[],[],[],opts);
        %[prs{p},fval(p)] = ga(@(x) VRdist_IF(x,stim,spvec,tc),...
        %    length(pp_all{p}),[],[],[],[],[],[],[],opts);
    end
end

[fval_min,idx] = min(fval);
fprintf(['Overall minimum: ' num2str(fval_min) ' (run ' num2str(idx) '); old minimum: ' num2str(fval_old) '\n'])

params_new = prs{idx};
