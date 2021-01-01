load('C:\Aneesha\StepIndent\afferent.mat')
load('C:\Users\somlab\Documents\Downloads\Nat_stim_filt.mat')
trace_now=Data.coffee_cup(:,1);
trace_now=trace_now(1200:4400);
%trace_now=trace_now([1:1136,2677:end]);
s = Stimulus(trace_now,[0 0],64,10);
r=a.response(s);

binwidth=.001;
for num_aff=1:length(r.responses)
    
    [n,bin]=histc(r.responses(num_aff).spikes,(-2:binwidth:100));
    spike_mat(num_aff,:)=n;
end
hold on


%Spike_Mat{PINS}{ampp}=spike_mat;se

%subplot(2,2,[1 3])
xlabel('Time (s)')
ylabel('Aggregate Afferent Spike Count (normalized)')
%xlim([-.15 1.15])

spike_matmax=spike_mat;
maxxFR= 1;        %max(sum(spike_matmax));
yyaxis left 
plot(-2:binwidth:100,smooth(sum(spike_mat))./maxxFR,'LineWidth',2);
xlabel('Time (s)')

hold on
timeindx=[-.03:1/(64*4):50.2];
timeindx=timeindx(1:length(trace_now));
yyaxis right
plot(timeindx,trace_now./max(trace_now));
xlim([0 51])
legend('aggregate firing', 'force trace')
legend boxoff

xlim([0 50])
% subplot(2,2,2)
% 
% xlabel('Time (s)')
% ylabel('Aggregate Afferent Spike Count (normalized)')
% %xlim([-.15 1.15])
% 
% spike_matmax=spike_mat;
% maxxFR=max(sum(spike_matmax));
% plot(-2:binwidth:100,sum(spike_mat)./maxxFR,'LineWidth',2);
% xlabel('Time (s)')
% 
% hold on
% plot([-.05:1/64:49.95],trace_now./max(trace_now));
% xlim([0 51])
% legend('aggregate firing', 'force trace')
% legend boxoff
% xlim([3 5])
% 
% 
% subplot(2,2,4)
% 
% xlabel('Time (s)')
% ylabel('Aggregate Afferent Spike Count (normalized)')
% %xlim([-.15 1.15])
% 
% spike_matmax=spike_mat;
% maxxFR=max(sum(spike_matmax));
% plot(-2:binwidth:100,sum(spike_mat)./maxxFR,'LineWidth',2);
% xlabel('Time (s)')
% 
% hold on
% plot([-.05:1/64:49.95],trace_now./max(trace_now));
% xlim([0 51])
% legend('aggregate firing', 'force trace')
% legend boxoff
% xlim([45 50])
