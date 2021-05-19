function ramp_and_hold()

close all;

%% show SA responses at different depths

amp = [0.3 0.54 0.78 1.02 1.44];
a = affpop_single_models();

figure(1)
set(1,'Position',[0 0 500 1000]);
for i=1:length(amp)
    subplot(length(amp),1,i)
    s = stim_ramp(amp(i),0.773,[0 0],5000,0.01,'lin',0.5);
    r(i) = a.response(s);
    for ii=1:sum(a.iSA1)
        plot_spikes(r(i).responses(findk(a.iSA1,ii)).spikes,'neuron_offset',ii,'hold','on')
    end
    ylim([1 5])
    xlim([-0.1 0.873])
end

screenshot(1,'figs/ramp_and_hold01_curr','pdf')
close(1)

%% show linearity of SA responses

amp = linspace(0,1.5,16);
for i=1:length(amp)
    s = stim_ramp(amp(i),0.6,[0 0],5000,0.01,'lin',0.5);
    r(i) = a.response(s);
    rr(:,i) = r(i).rate(a.iSA1);
end

figure(2)
plot(amp,rr','o')
hold on
lsline

screenshot(2,'figs/ramp_and_hold02_curr','pdf')
close(2)

%% onset and offset transients for RAs

ramp_vel = [60 26 15 9 6 3 2.2 1.7 1.4 0.95 0.78 0.55];
ramp_len = 0.5./ramp_vel;

figure(3)
set(3,'Position',[0 0 900 1000]);
idx = reshape(1:12,3,[])';
idx = idx(:);
for i=1:length(ramp_vel)
    subplot(4,3,idx(i))
    hold on
    s = stim_ramp(0.5,2*max(ramp_len),[0 0],5000,ramp_len(i),'lin',1);
    r(i) = a.response(s);
    for ii=1:sum(a.iRA)
        plot_spikes(r(i).responses(findk(a.iRA,ii)).spikes,'neuron_offset',ii,'hold','on')
    end
    ylim([1 10])
    xlim([-0.1 2*max(ramp_len)+0.1])
    title([num2str(ramp_vel(i)) 'mm/s'])
end

screenshot(3,'figs/ramp_and_hold03_curr','pdf')
close(3)
