function h = hist_spikes(spikes,edges)

h = histc(spikes,edges,1);
if isempty(h)
    h = zeros(length(edges),1);
end
h(end) = [];
