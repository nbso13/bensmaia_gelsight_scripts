function [idx,aff_class] = find_common_neurons(varargin)

rej = {'pma_15_03','pma_13_00','pma_09_03','pma_13_01','pma_13_02','pma_14_00'};

nid = varargin{1}.neuron_id;
aff_class = varargin{1}.nType;
idx = NaN*zeros(length(varargin),length(nid));
idx(1,:) = 1:length(nid);

for i=1:length(nid)
    for s=2:length(varargin)
        pos = strmatch(nid{i},varargin{s}.neuron_id);
        if ~isempty(pos)
            idx(s,i) = pos;
        end
    end
end

% remove rejected neurons
for i=1:length(rej)
    ii = strmatch(rej{i},nid);
    idx(:,ii) = NaN;
end
