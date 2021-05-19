% adaptation.m

a = affpop_single_models();

amp = linspace(0,1.5,16);
for i=1:length(amp)
    s = stim_ramp(amp(i),0.6,[0 0],5000,0.01,'lin',0.5);
    r(i) = a.response(s);
    rr(:,i) = r(i).rate(a.iSA1);
end

figure(1)
set(gcf,'pos',[0 0 240 474]);
subplot(211)

plot(amp,rr','o')
hold on
lsline
box off
ylim([0 Inf])

xlabel('Indentation depth [mm]')
ylabel('Firing rate [Hz]')

screenshot(1,'adaptation1','pdf')
%close(1)
