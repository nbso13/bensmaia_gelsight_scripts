% model_fits.m
shouldUseGLM = 0;

a = affpop_single_models();
if(exist('shouldUseGLM','var') && shouldUseGLM)
    a = affpop_single_models([],[],'Model','GLM');
end

%% Sinusoids

load SINE.mat

%      SA         RA                         PC
idx = [2 6 7 10   1 3 13 14 15 16 22 23 25   8 9 11 27];
num = length(idx);

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
    for n=1:num
        r_act(n,i) = length(SINE.ObservedSpikes{idx(n)}{idx_sine(i),1})/...
            stim_len(i);
        times_act(n,i)=SINE.ObservedSpikes{idx(n)}(idx_sine(i),1);
    end
    
    s = Stimulus(desiredsinepos{i}'/1000,[0 0],20000,0.5);
    r = a.response(s);
    
    r_sim(:,i) = r.rate;
    times_sim(:,i)={r.responses(:).spikes};
end

% compute timing specific variables
% number of spikes
n_act=cellfun(@length,times_act);
n_sim=cellfun(@length,times_sim);
nidx=[1:4 1:9 1:4];

% ISI
isi_act=cellfun(@diff,times_act,'uni',0);
isi_sim=cellfun(@diff,times_sim,'uni',0);

% median ISI in a trial
medisi_act=cellfun(@median,isi_act);
medisi_sim=cellfun(@median,isi_sim);

nn_act=~isnan(medisi_act); % non-nan indices
nn_sim=~isnan(medisi_sim);

% thresholds
for f=1:length(freqs_unique)
    idx = find(freqs==freqs_unique(f));
    medisi_act_f = medisi_act(:,idx);
    medisi_sim_f = medisi_sim(:,idx);
    for n=1:size(medisi_act,1)
        min_act_f = min(amps(idx(~isnan(medisi_act_f(n,:)))));
        if isempty(min_act_f)
            thres_act(n,f) = NaN;
        else
            thres_act(n,f) = min_act_f;
        end
        min_sim_f = min(amps(idx(~isnan(medisi_sim_f(n,:)))));
        if isempty(min_sim_f)
            thres_sim(n,f) = NaN;
        else
            thres_sim(n,f) = min_sim_f;
        end
    end
end

%% Noise

load NOISE.mat

%      SA        RA                         PC
idx = [2 4 5 8   1 3 11 12 13 14 18 19 20   6 7 9 21];

try
    load desirednoisepos.mat
catch
    process_noise_pos;
    load desirednoisepos.mat
end

for i=1:60
    for n=1:num
        r_act_noise(n,i) = length(NOISE.ObservedSpikes{idx(n)}{i,1});
    end
    
    s = Stimulus(reshape(desirednoisepos{i},[],1)/1000,[0 0],20000,0.5);
    r = a.response(s);
    r_resp = r.responses;
    
    r_sim_noise(:,i) = r.rate;
end

%% Diharms

load DI.mat

%      SA         RA                         PC
idx = [2 6 7 10   1 3 13 14 15 16 21 22 24   8 9 11 27];
num = length(idx);

try
    load desireddiharmpos.mat
catch
    process_diharm_pos;
    load desireddiharmpos.mat
end

for i=1:size(stim,1)
    stim_ind = find(DI.state_v(:,5)==stim(i,1) & DI.state_v(:,6)==stim(i,2)...
        & DI.state_v(:,2)==stim(i,3));
    for n=1:num
        r_act_di(n,i) = length(DI.ObservedSpikes{idx(n)}{stim_ind,1})/...
            DI.state_v(stim_ind,end);
    end
    
    s = Stimulus(desireddiharmpos{i}/1000,[0 0],20000,0.5);
    r = a.response(s);
    r_sim_di(:,i) = r.rate;
end

%% Actual vs predicted rates: TRAINING DATA

figure(1)
set(1,'pos',[1 1 1420 1080],'PaperPosition', [.25 .25 25 25]);
s_idx = [1 2 3 4 6 7 8 9 10 11 12 13 14 16 17 18 19];
for n=1:num
    subplot(4,5,s_idx(n))
    hold on
    if strcmpi(a.afferents(n).class,'PC')
        F = 251;
    else
        F = 120;
    end
    cc = corrcoef(r_act(n,1:F),r_sim(n,1:F),'rows','complete');
    r2_old(n) = 1 - sum((r_act(n,1:F)-r_sim(n,1:F)).^2)/sum((r_act(n,1:F)-mean(r_act(n,1:F))).^2);
    
    stats(n) = regstats(r_act(n,1:F),r_sim(n,1:F),'linear');
    r2(n) = stats(n).rsquare;
    intercept(n) = stats(n).beta(1);
    slope(n) = stats(n).beta(2);
    [p_intercept(n),F_intercept(n)] = linhyptest(stats(n).beta,stats(n).covb,0,[1 0],stats(n).fstat.dfe);
    [p_slope(n),F_slope(n)] = linhyptest(stats(n).beta,stats(n).covb,1,[0 1],stats(n).fstat.dfe);
    
    sse(n) = sum((r_act(n,1:F)-r_sim(n,1:F)).^2);
    Ubias(n) = (F*(mean(r_act(n,1:F))-mean(r_sim(n,1:F))).^2)/sse(n);
    Uslope(n) = (((slope(n)-1).^2)*sum((r_sim(n,1:F)-mean(r_sim(n,1:F))).^2))/sse(n);
    Uerror(n) = sum((stats(n).yhat'-r_act(n,1:F)).^2)/sse(n);
    
    plot(r_sim(n,1:F),r_act(n,1:F),'o','MarkerSize',5,'Color',[0 0.4470 0.7410]);
    
    xx = get(gca,'xlim');
    plot(xx,intercept(n)+slope(n)*xx,'Color',[0 0.4470 0.7410])
    
    if strcmpi(a.afferents(n).class,'PC')
        F = 1;
    else
        F = 41;
    end
    cc_noise = corrcoef(r_act_noise(n,F:end),r_sim_noise(n,F:end),'rows','complete');
    r2_noise_old(n) = 1 - sum((r_act_noise(n,F:end)-r_sim_noise(n,F:end)).^2)/sum((r_act_noise(n,F:end)-mean(r_act_noise(n,F:end))).^2);
    
    stats_noise(n) = regstats(r_act_noise(n,F:end),r_sim_noise(n,F:end),'linear');
    r2_noise(n) = stats_noise(n).rsquare;
    intercept_noise(n) = stats_noise(n).beta(1);
    slope_noise(n) = stats_noise(n).beta(2);
    [p_intercept_noise(n),F_intercept_noise(n)] = linhyptest(stats_noise(n).beta,stats_noise(n).covb,0,[1 0],stats_noise(n).fstat.dfe);
    [p_slope_noise(n),F_slope_noise(n)] = linhyptest(stats_noise(n).beta,stats_noise(n).covb,1,[0 1],stats_noise(n).fstat.dfe);
    
    sse_noise(n) = sum((r_act_noise(n,F:end)-r_sim_noise(n,F:end)).^2);
    Ubias_noise(n) = (length(r_act_noise(n,F:end))*(mean(r_act_noise(n,F:end))-mean(r_sim_noise(n,F:end))).^2)/sse_noise(n);
    Uslope_noise(n) = (((slope_noise(n)-1).^2)*sum((r_sim_noise(n,F:end)-mean(r_sim_noise(n,F:end))).^2))/sse_noise(n);
    Uerror_noise(n) = sum((stats_noise(n).yhat'-r_act_noise(n,F:end)).^2)/sse_noise(n);
    
    plot(r_sim_noise(n,F:end),r_act_noise(n,F:end),'o','MarkerSize',5,'Color',[0.8500 0.3250 0.0980])
    plot(xx,intercept_noise(n)+slope_noise(n)*xx,'Color',[0.8500 0.3250 0.0980])
    
    unityslope;
    text(0.1,0.95,['R^2 = ' num2str(r2(n),2)],'units','normalized','Color',[0 0.4470 0.7410])
    text(0.1,0.85,['R^2 = ' num2str(r2_noise(n),2)],'units','normalized','Color',[0.8500 0.3250 0.0980])
    box off
    
    xlim([0 Inf])
    ylim([0 Inf])
    
    axis square
    
    %     statsX(n) = regstats([r_act(n,1:F) r_act_noise(n,F:end)],[r_sim(n,1:F) r_sim_noise(n,F:end)],'linear');
    %     r2X(n) = statsX(n).rsquare;
    %     interceptX(n) = statsX(n).beta(1);
    %     slopeX(n) = statsX(n).beta(2);
    %     [p_interceptX(n),F_interceptX(n)] = linhyptest(statsX(n).beta,statsX(n).covb,0,[1 0],statsX(n).fstat.dfe);
    %     [p_slopeX(n),F_slopeX(n)] = linhyptest(statsX(n).beta,statsX(n).covb,1,[0 1],statsX(n).fstat.dfe);
end

subplot(4,5,16)
xlabel('Predicted rate [Hz]')
ylabel('Measured rate [Hz]')
legend('Sine','Noise','Location','SouthEast')

screenshot(1,'model_fits1','pdf')
if(screenshot)
    close(1)
end

%% Actual vs predicted rates: TEST DATA

figure(2)
set(2,'pos',[1 1 1420 1080],'PaperPosition', [.25 .25 25 25]);
s_idx = [1 2 3 4 6 7 8 9 10 11 12 13 14 16 17 18 19];
for n=1:num
    subplot(4,5,s_idx(n))
    hold on
    
    if strcmpi(a.afferents(n).class,'PC')
        F = [1:138];
    else
        F = [1:55 76:80];
    end
    cc_di = corrcoef(r_act_di(n,F),r_sim_di(n,F));
    r2_di_old(n) = 1 - sum((r_act_di(n,F)-r_sim_di(n,F)).^2)/sum((r_act_di(n,F)-mean(r_act_di(n,F))).^2);
    
    stats_di(n) = regstats(r_act_di(n,F),r_sim_di(n,F),'linear');
    r2_di(n) = stats_di(n).rsquare;
    intercept_di(n) = stats_di(n).beta(1);
    slope_di(n) = stats_di(n).beta(2);
    [p_intercept_di(n),F_intercept_di(n)] = linhyptest(stats_di(n).beta,stats_di(n).covb,0,[1 0],stats_di(n).fstat.dfe);
    [p_slope_di(n),F_slope_di(n)] = linhyptest(stats_di(n).beta,stats_di(n).covb,1,[0 1],stats_di(n).fstat.dfe);
    
    sse_di(n) = sum((r_act_di(n,F)-r_sim_di(n,F)).^2);
    Ubias_di(n) = (length(r_act_di(n,F))*(mean(r_act_di(n,F))-mean(r_sim_di(n,F))).^2)/sse_di(n);
    Uslope_di(n) = (((slope_di(n)-1).^2)*sum((r_sim_di(n,F)-mean(r_sim_di(n,F))).^2))/sse_di(n);
    Uerror_di(n) = sum((stats_di(n).yhat'-r_act_di(n,F)).^2)/sse_di(n);
    
    plot(r_sim_di(n,F),r_act_di(n,F),'o','MarkerSize',5,'Color',[0.4940 0.1840 0.5560])
    
    xx = get(gca,'xlim');
    plot(xx,intercept_di(n)+slope_di(n)*xx,'Color',[0.4940 0.1840 0.5560])
    
    unityslope;
    text(0.1,0.95,['R^2 = ' num2str(r2_di(n),2)],'units','normalized','Color',[0.4940 0.1840 0.5560])
    box off
    
    xlim([0 Inf])
    ylim([0 Inf])
    
    axis square
end

subplot(4,5,16)
xlabel('Predicted rate [Hz]')
ylabel('Measured rate [Hz]')
legend('Diharm','Location','SouthEast')

screenshot(2,'model_fits2','pdf')
if(screenshot)
    close(2)
end
%% export all regression values

T = table();

% SINE
T.R2_sine = r2';

T.slope_sine = num2str(slope','%1.2f');
T.slope_sine = [T.slope_sine repmat(char(32),17,2)];
T.slope_sine(p_slope<0.05,end-1) = '*';
T.slope_sine(p_slope<0.01,end) = '*';

T.intercept_sine = num2str(intercept','%1.2f');
T.intercept_sine = [T.intercept_sine repmat(char(32),17,2)];
T.intercept_sine(p_intercept<0.05,end-1) = '*';
T.intercept_sine(p_intercept<0.01,end) = '*';

T.Ubias_sine = Ubias';
T.Uconsistency_sine = Uslope';
T.Uerror_sine = Uerror';

% NOISE
T.R2_noise = r2_noise';

T.slope_noise = num2str(slope_noise','%1.2f');
T.slope_noise = [T.slope_noise repmat(char(32),17,2)];
T.slope_noise(p_slope_noise<0.05,end-1) = '*';
T.slope_noise(p_slope_noise<0.01,end) = '*';

T.intercept_noise = num2str(intercept_noise','%1.2f');
T.intercept_noise = [T.intercept_noise repmat(char(32),17,2)];
T.intercept_noise(p_intercept_noise<0.05,end-1) = '*';
T.intercept_noise(p_intercept_noise<0.01,end) = '*';

T.Ubias_noise = Ubias_noise';
T.Uconsistency_noise = Uslope_noise';
T.Uerror_noise = Uerror_noise';

% DIHARM
T.R2_diharm = r2_di';

T.slope_diharm = num2str(slope_di','%1.2f');
T.slope_diharm = [T.slope_diharm repmat(char(32),17,2)];
T.slope_diharm(p_slope_di<0.05,end-1) = '*';
T.slope_diharm(p_slope_di<0.01,end) = '*';

T.intercept_diharm = num2str(intercept_di','%1.2f');
T.intercept_diharm = [T.intercept_diharm repmat(char(32),17,2)];
T.intercept_diharm(p_intercept_di<0.05,end-1) = '*';
T.intercept_diharm(p_intercept_di<0.01,end) = '*';

T.Ubias_diharm = Ubias_di';
T.Uconsistency_diharm = Uslope_di';
T.Uerror_diharm = Uerror_di';

% write results to file
writetable(T,'model_fits.xlsx')

return
%% figure for talk
R_act{1,1}=r_act(1:4,1:120);
R_act{2,1}=r_act(5:13,1:120);
R_act{3,1}=r_act(14:17,:);

R_act{1,2}=r_act_noise(1:4,41:end);
R_act{2,2}=r_act_noise(5:13,41:end);
R_act{3,2}=r_act_noise(14:17,:);

R_act{1,3}=r_act_di(1:4,[1:55 76:80]);
R_act{2,3}=r_act_di(5:13,[1:55 76:80]);
R_act{3,3}=r_act_di(14:17,1:138);

R_sim{1,1}=r_sim(1:4,1:120);
R_sim{2,1}=r_sim(5:13,1:120);
R_sim{3,1}=r_sim(14:17,:);

R_sim{1,2}=r_sim_noise(1:4,41:end);
R_sim{2,2}=r_sim_noise(5:13,41:end);
R_sim{3,2}=r_sim_noise(14:17,:);

R_sim{1,3}=r_sim_di(1:4,[1:55 76:80]);
R_sim{2,3}=r_sim_di(5:13,[1:55 76:80]);
R_sim{3,3}=r_sim_di(14:17,1:138);

dm=@(x) bsxfun(@times,mean(x)',[1 1])';
pmsd=@(x) bsxfun(@plus,mean(x)',bsxfun(@times,sem(x)'/2,[-1 1]))';
rsquare=@(actu,pred) 1-nansum((actu-pred).^2)/nansum((actu-mean(actu)).^2);
path='D:\Dropbox (INMACOSY)\confs\2017_06_WH\';

type={'SA1','RA','PC'};
stim={'sine','noise','diharm'};
lims=[80 150 680];
cols=[0,0,0;.9,.1,.1;.1,.3,.9];
R=zeros(3,3);
fh=newfig(['rates ' type{ii}],'pos',[ii*400 350 270*3 270]/40);
pl=[.24 .52 .80]

for jj=1:3 % stimulus
    for ii=1:3 % class
        ax(ii)=subplot(1,3,ii);
        set(ax(ii),'nextplot','add')
        plot(ax(ii),[0 lims(ii)],[0 lims(ii)],'-','col',[.6 .6 .6],'linew',1.5)
        
        R(ii,jj)=rsquare(mean(R_act{ii,jj}),mean(R_sim{ii,jj}));
        options={'linew',1,'col',cols(jj,:),'markers',3.5};
        plot(ax(ii),pmsd(R_sim{ii,jj}),dm(R_act{ii,jj}),options{:});
        plot(ax(ii),dm(R_sim{ii,jj}),pmsd(R_act{ii,jj}),options{:});
        h(ii,jj)=plot(ax(ii),mean(R_sim{ii,jj}),mean(R_act{ii,jj}),'s',...
            options{:},'displayname',[stim{jj} ' (R^2=' num2str(R(ii,jj),2) ')']);
        set(h(ii,jj),'markeredgec','none','markerfacec',cols(jj,:))
        
        
        title(ax(ii),type{ii})
        hl=lgd(ax(ii),h(ii,1:jj),'location','SE'); set(hl,'pos',[pl(ii) .2 0.1855 0.2260])
        axis(ax(ii),'equal')
        xt=get(ax(ii),'xtick');
        set(ax(ii),'xlim',[0 lims(ii)],'ylim',[0 lims(ii)],'ytick',xt)
        
        ylabel(ax(ii),'Measured rate [Hz]')
        xlabel(ax(ii),'Predicted rate [Hz]')
        
    end
    if(jj==2)
        screenshot(fh,[path 'talk_fitrate' type{ii} 'without_diharm'])
    elseif(jj==3)
        screenshot(fh,[path 'talk_fitrate' type{ii} 'with_diharm'])
    end
end



