% example_fits.m

idx = [1 9 4];

cols=affcol;

%%
a = AfferentPopulation();
a.add_afferents('SA1',[0 0],'idx',idx(1));
a.add_afferents('RA',[0 0],'idx',idx(2));
a.add_afferents('PC',[0 0],'idx',idx(3));

load SINE.mat
stim_len = SINE.state_v(:,8);
idx_sine = 1:length(stim_len);

[freqs,ind_sorted] = sort(SINE.state_v(idx_sine,5),'ascend');
amps=SINE.state_v(idx_sine(ind_sorted),2);
ntypes=SINE.nType(idx);
freqs_unique = unique(freqs);
idx_sine = idx_sine(ind_sorted);
stim_len = stim_len(ind_sorted);

try
    load desiredsinepos.mat
catch
    process_sine_pos;
    load desiredsinepos.mat
end

for i=1:length(idx_sine)
    
    s = Stimulus(desiredsinepos{i}'/1000,[0 0],20000,0.5);
    r = a.response(s);
    
    r_sim(:,i) = r.rate;
    times_sim(:,i)={r.responses(:).spikes};
    
end

% compute timing specific variables
% number of spikes
n_sim=cellfun(@length,times_sim);

% ISI
isi_sim=cellfun(@diff,times_sim,'uni',0);

% median ISI in a trial
medisi_sim=cellfun(@median,isi_sim);

nn_sim=~isnan(medisi_sim);

% thresholds
for f=1:length(freqs_unique)
    idx = find(freqs==freqs_unique(f));
    medisi_sim_f = medisi_sim(:,idx);
    for n=1:size(medisi_sim,1)
        min_sim_f = min(amps(idx(~isnan(medisi_sim_f(n,:)))));
        if isempty(min_sim_f)
            thres_sim(n,f) = NaN;
        else
            thres_sim(n,f) = min_sim_f;
        end
    end
end


nneur=[3];
nnumb=[ones(1,3)];
nidxs=[1:3];

% plot all
figure(1)
set(1,'pos',[0 0 240 711])


colormap(1,flipud(jet))
tri=delaunay(freqs,amps);

for id=1:3
    subplot(3,1,id)
    hold on
    
    f=nnumb(id);
    n=nidxs(id);
    %set(ax{f,n},'nextplot','add');
    % plot empty circles
    scatter(freqs,amps,40,[.7 .7 .7])
    % plot coloredfilled dots
    idx2=find(~isnan(medisi_sim(id,:)));
    scatter(freqs(idx2),amps(idx2),40,log10(medisi_sim(id,idx2)),'filled')
    % draw a line for thresholds
    plot(freqs_unique,thres_sim(id,:),'k','LineWidth',1.5)
    
    
    set(gca,'xscale','log','yscale','log','xlim',[1 1000],...
        'ylim',[.0625 500],'clim',[-3 0],...
        'xtick',[1 10 100 1000],'ytick',[1 10 100],'xminorgrid','off',...
        'yminorgrid','off');
    xlabel('Frequency [Hz]')
    ylabel('Amplitude [um]')
end

%h(1)=colorbar();
%set(gca,'visible','off','clim',[-3 0])
%set(h,'ticks',[-3 -2 -1 0],'ticklabels',[.001 .01 .1 1]);
%h(1).Label.String='ISI [s]';

screenshot(1,'example_fits1','pdf')
close(1)

%% Spike trains for ramps and sinusoids

figure(2)
set(2,'pos',[0 0 240 900])

s(1) = stim_ramp(0.5,0.5,[0 0],5000,0.01,'sine',1);
s(2) = stim_ramp(0.75,0.5,[0 0],5000,0.075,'sine',1);
s(3) = stim_sine(200,0.005,[],0.15,[],[],0.02);

num_plots = length(s)*4;

c = 1;
for ss =1:length(s)
    r = a.response(s(ss));
    subplot(num_plots,1,c)
    plot(linspace(0,s(ss).duration,length(s(ss).trace)),s(ss).trace,'k','LineWidth',1.5)
    xlim([0 s(ss).duration])
    if c==1 || c==5
        ylim([0 1])
    end
    box off
    c = c + 1;
    for id=1:3
        subplot(num_plots,1,c)
        hold on
        plot_spikes(r.responses(id).spikes,'linewidth',1.5,'color',cols(id,:))
        xlim([0 s(ss).duration])
        c = c + 1;
    end
end

screenshot(2,'example_fits2','pdf')
close(2)

%% plot receptive field examples

rads = [1.2754 1.9287 9.0820;
        5.3037 4.2725 41.9678];

figure(3)
[origin,~,pxl_per_mm,regionprop] = plot_hand('axes',false,'names',false,'centers',false);
clf

c = 1;
for i=1:2
    for j=1:3
        subplot(2,3,c)
        hold on
        plot(regionprop(3).Boundary(:,1),regionprop(3).Boundary(:,2),'k','LineWidth',1.5)
        plot(regionprop(5).Boundary(:,1),regionprop(5).Boundary(:,2),'k','LineWidth',1.5)
        viscircles(origin,rads(i,j),'Color',cols(j,:));
        xlim([100 200])
        ylim([250 500])
        axis equal
        axis off
        
        c = c + 1;
    end
end

screenshot(3,'example_fits3','pdf')
close(3)
