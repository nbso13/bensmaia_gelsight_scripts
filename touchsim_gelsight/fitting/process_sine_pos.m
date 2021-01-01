function process_sine_pos

try
    load SINE.mat
catch
    try
        load \\bsl-somsrv1\data\processed\per\pma\L2_str\SINE.mat
    catch
        error('Unable to find SINE.mat locally or on server.')
    end
end

[all_freqs,ind_sorted] = sort(SINE.state_v(:,5),'ascend');
all_amps = SINE.state_v(ind_sorted,2);
all_stim_len = SINE.state_v(ind_sorted,8);

freqs = unique(all_freqs);

counter = 1;
for f=1:length(freqs)
    trace = [];
    for i=1:15
        file = fopen(['sine\sine_' num2str(freqs(f)) '_' num2str(i) '.bin']);
        if file==-1
            file = fopen(['\\bsl-somsrv1\data\processed\per\pma\desired\sine\sine_' num2str(freqs(f)) '_' num2str(i) '.bin']);
        end
        if file==-1
            error('Unable to find sineXX.bin files locally or on server.')
        end
        
        trace(i,:) = fread(file,Inf,'double');
        fclose(file);
    end
    des_amps = max(trace,[],2);
    amps = all_amps(all_freqs==freqs(f));
    for a=1:length(amps)
        [d,ind] = sort((des_amps-amps(a)).^2,'ascend');
        diff(counter) = d(1);
        desiredsinepos{counter} = trace(ind(1),1002:end-1002);
        counter = counter + 1;
    end
end

save desiredsinepos.mat desiredsinepos
