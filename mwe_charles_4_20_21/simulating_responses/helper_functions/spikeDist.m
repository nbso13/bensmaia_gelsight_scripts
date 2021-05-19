
function dist = spikeDist(spike1, spike2, q)
ignore = 0;

if xor(length(spike1) == 0,  length(spike2) == 0)
    if length(spike1) == 0
        if isnan(spike2)
            ignore = 1;
        end
    else
        if isnan(spike1)
            ignore = 1;
        end
    end
end

if ignore == 1
    dist = nan;
else
    dist = spkdl([spike1; spike2],[1 length(spike1)+1], [length(spike1) length(spike1)+length(spike2)],q);
    dist = dist(2,:); % THIS LINE ADDED BY KHL
end

% the above if statements correct a problem that arose anytime one input
% was an empty vector and the other a nan. this combination should now
% always result in a nan output instead of a numerical distance.