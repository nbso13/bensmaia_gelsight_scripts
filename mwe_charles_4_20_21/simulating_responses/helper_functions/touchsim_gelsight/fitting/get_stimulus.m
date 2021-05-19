function [stim,num_stim] = get_stimulus(pad,varargin)
% [stim,num_stim] = get_stimulus(pad,varargin)

num_stim = 0;

for s=1:length(varargin)
    stimi{s} = varargin{s};
    num_stimi = size(stimi{s},1);
    num_stim = num_stim + num_stimi;
    
    stimi{s} = [stimi{s} zeros(num_stimi,pad)]';
    stimi{s} = reshape(stimi{s},[],1);    
end

stim = cell2mat(stimi');
%stim = smooth(stim,30);
