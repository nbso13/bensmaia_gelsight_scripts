clear all;close all;clc


ap=affpop_hand('D2d');

%% A letter
h = figure('Visible','off','pos',[100 100 300 300]);
axis off
text(0.5,0.5,'A','FontSize',250,'horiz','center')

A = print2array(h);
close(h)
A = mean(A,3);
A = A-min(A(:));
A = A/max(A(:));

wid = 10.1;
pins_per_mm=10;

[a,b] = size(A);
[x,y] = meshgrid(linspace(-wid/2/a*b,wid/2/a*b,b),linspace(-wid/2,wid/2,a));
[X,Y] = meshgrid(linspace(-wid/2/a*b,wid/2/a*b,ceil(wid/a*b*pins_per_mm)),linspace(-wid/2,wid/2,ceil(wid*pins_per_mm)));
A = interp2(x,y,A',X,Y);
A=1-A;

%%
sr=stim_ramp;
trace=bsxfun(@times,sr.trace,A(:)');
stim=Stimulus(trace,[X(:) Y(:)],sr.sampling_frequency);
rc=ap.response(stim);
spikes={rc.responses(:).spikes};
rates=cellfun(@length,spikes);

%%
% plot_spikes(spikes)

%%
subplot(131)
scatter(ap.location(ap.iSA1,2),-ap.location(ap.iSA1,1),[],rates(ap.iSA1),'filled')
subplot(132)
scatter(ap.location(ap.iRA,2),-ap.location(ap.iRA,1),[],rates(ap.iRA),'filled')
subplot(133)
scatter(ap.location(ap.iPC,2),-ap.location(ap.iPC,1),[],rates(ap.iPC),'filled')
colormap hot