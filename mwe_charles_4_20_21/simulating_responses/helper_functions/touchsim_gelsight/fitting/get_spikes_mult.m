function [sp,spvec] = get_spikes_mult(L2_str,aff_ind,idx,params,rep)
% [sp,spvec] = get_spikes_mult(L2_str,aff_ind,idx,params,rep)

if nargin<5
    rep = 1;
end

for s=1:length(L2_str)
        
    edges = cellfun(@linspace,mat2cell(params{s}(:,1),ones(length(idx{s}),1)),...
        mat2cell(params{s}(:,1)+params{s}(:,2),ones(length(idx{s}),1)),...
        mat2cell(round(params{s}(:,2)*5000)+1,ones(length(idx{s}),1)),'UniformOutput',false);
    hspikes = cellfun(@hist_spikes,L2_str{s}.ObservedSpikes{aff_ind{s}}(idx{s},rep),...
        edges,'UniformOutput',false);
    hspikes_con = cellfun(@vertcat,hspikes,...
        mat2cell(zeros(900,length(idx{s})),900,ones(1,length(idx{s})))','UniformOutput',false);
    
    spvec{s} = cell2mat(hspikes_con);
end
spvec = cell2mat(spvec');
sp = find(spvec);
