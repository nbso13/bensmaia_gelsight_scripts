function process_diharm_pos

try
    load DI.mat
catch
    try
        load \\bsl-somsrv1\data\processed\per\pma\L2_str\DI.mat
    catch
        error('Unable to find DI.mat locally or on server.')
    end
end

[all_freqs,ind_sorted] = sortrows(DI.state_v(:,5:6));
all_amps = DI.state_v(ind_sorted,2);
all_stim_len = DI.state_v(ind_sorted,8);

freqs = unique(all_freqs,'rows');

counter = 1;
for f=1:size(freqs,1)
    amps_sub = sort(all_amps(all_freqs(:,1)==freqs(f,1) & all_freqs(:,2)==freqs(f,2)),'ascend');
    trace = [];
    for i=1:length(amps_sub)
        try
            error
            file = fopen(['diharm\diharm_' num2str(freqs(f,1)) '_' num2str(freqs(f,2)) '_' num2str(i) '_0.bin']);
        catch
            try
                file = fopen(['\\bsl-somsrv1\data\processed\per\pma\desired\diharm\diharm_' num2str(freqs(f,1)) '_' num2str(freqs(f,2)) '_' num2str(i) '_0.bin']);
            catch
                error('Unable to find diharmXX.bin files locally or on server.')
            end
        end
        
        trace = fread(file,Inf,'double');
        desireddiharmpos{counter} = trace(1001:end-1001);
        stim(counter,1) = freqs(f,1);
        stim(counter,2) = freqs(f,2);
        stim(counter,3) = amps_sub(i);
        fclose(file);
        counter = counter + 1;
    end
end

save desireddiharmpos.mat desireddiharmpos stim
