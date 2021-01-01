function plot_fit(stim_expanded,spvec,predvec,vmem)
% plot_fit(stim_expanded,spvec,predvec)

%xx = linspace(0,size(stim_expanded,1)/5,size(stim_expanded,1));
xx = linspace(0,1,length(spvec));

clf
hold on
%plot(xx,stim_expanded(:,1))
%plot(xx,-100*stim_expanded(:,5))
plot(xx,20*spvec,'r')
plot(xx,-20*predvec,'g')
if nargin==4
    plot(xx,vmem,'k')
end
drawnow
