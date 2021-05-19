% setup and loading data
ccc
datasets={'sine','noise','diharm'};

% RELOAD OR COMPUTE REPONSES (takes some time)
flag_reload=1;

% load measured spiking data
[data{1},data{2},data{3}]=getVibData;

ntnum=[ones(4,1);ones(9,1)*2;ones(4,1)*3];
aff={'SA1','RA','PC'};
datestr=char(datetime('now','TimeZone','local','Format','dd_MM_yyyy_HH_mm_ss'));

% SETUP: 1 for sine, 2 for noise, 3 for diharm
for ds=3
    fprintf('Dataset: %s \n',datasets{ds})
    fn=fieldnames(data{ds});
    for ii=1:length(fn), eval([fn{ii} '=data{ds}.' fn{ii} ';']); end
    switch ds % pos files
        case 1, posdata=load('mat/desiredsinepos.mat');
            pos=posdata.desiredsinepos; clear posdata; toff=0.05;
        case 2, posdata=load('mat/desirednoisepos.mat');
            pos=posdata.desirednoisepos; clear posdata; toff=0.5;
        case 3, posdata=load('mat/desireddiharmpos.mat');
            pos=posdata.desireddiharmpos; clear posdata;  toff=0.05;
    end
    for ii=1:naff
        spikes{ii}=cellfun(@(x) x-toff,spikes{ii},'uni',0);
        spikes{ii}=cellfun(@(x) x(x>0),spikes{ii},'uni',0);
    end
    spk_act=spikes;
    %% model responses
    if(~flag_reload) % reload model responses if already computed
        fprintf('computing model resp: '); tic
        a_IF = affpop_single_models([],[],'noisy',true);
        a_GLM = affpop_single_models([],[],'model','GLM');
        spk_IF=cell(naff,1);
        spk_GLM=cell(naff,1);
        for ii=1:ncond
            s = Stimulus(reshape(pos{ii},[],1)/1000,[0 0],20000,0.5);
            for jj=1:nrep
                r_IF = {a_IF.response(s).responses(:).spikes};
                r_GLM = {a_GLM.response(s).responses(:).spikes};
                for kk=1:naff
                    spk_IF{kk}{ii,jj}=r_IF{kk};
                    spk_GLM{kk}{ii,jj}=r_GLM{kk};
                end
            end
        end
        toc
        % save
        filename=[datasets{ds} '_' datestr];
        save(['mat/timing/' filename],'spk_IF','spk_GLM')
        fprintf('resp saved \n');
    else
        fprintf('reloading model resp: '); tic
        files=dir(['mat/timing/' datasets{ds} '*.mat']);
        filename=files(end).name;
        load(['mat/timing/' filename]);
        toc
    end
    
    %% computing ISI histograms
    fprintf('computing ISI histograms: '); tic
    isi_bin=logspace(-3.2,0,40);
    binwidth = diff(isi_bin);
    edges = [-Inf, isi_bin(1:end-1)+binwidth/2, Inf];
    isi_act=cell(naff,1);      isih_act=cell(ncond+1,naff);
    isi_IF=cell(naff,1);       isih_IF=cell(ncond+1,naff);
    isi_GLM=cell(naff,1);      isih_GLM=cell(ncond+1,naff);
    for ii=1:naff
        isi_act{ii}=cellfun(@diff,spk_act{ii},'uni',0);
        lisi_act=num2cell(cellfun(@(x) histcountsmex(x,edges),isi_act{ii},'uni',0),2);
        isih_act(1:end-1,ii)=cellfun(@(x) cat(1,x{:})',lisi_act,'uni',0);
        isih_act{end,ii}=sum(cat(3,isih_act{:,ii}),3);
        
        isi_IF{ii}=cellfun(@diff,spk_IF{ii},'uni',0);
        lisi_IF=num2cell(cellfun(@(x) histcountsmex(x,edges),isi_IF{ii},'uni',0),2);
        isih_IF(1:end-1,ii)=cellfun(@(x) cat(1,x{:})',lisi_IF,'uni',0);
        isih_IF{end,ii}=sum(cat(3,isih_IF{:,ii}),3);
        
        isi_GLM{ii}=cellfun(@diff,spk_GLM{ii},'uni',0);
        lisi_GLM=num2cell(cellfun(@(x) histcountsmex(x,edges),isi_GLM{ii},'uni',0),2);
        isih_GLM(1:end-1,ii)=cellfun(@(x) cat(1,x{:})',lisi_GLM,'uni',0);
        isih_GLM{end,ii}=sum(cat(3,isih_GLM{:,ii}),3);
    end
    toc
    %% computing spike distances
    if ~exist(['mat/timing/spkd_' filename],'file')
        fprintf('computing spike distances: '); tic
        % select subset of conditions models have been fitted on
        conds={1:120 1:120 1:251; ... sines
            41:60 41:60 1:60; ... noises
            [1:55 76:80] [1:55 76:80] 1:138}; % diharms
        % parameters
        q=logspace(0.1,3.9,50); nq=length(q);
        lags=(-.005:.00025:.005)';
        spknumlim=3;
        % compute
        dActInt=cell(naff,1);     dIF=cell(naff,1);     dGLM=cell(naff,1);
        for ii=1:naff
            ldActInt=nan(ncond,nq); ldIF=nan(ncond,nq);  ldGLM=nan(ncond,nq);
            lspk_act=spk_act{ii};   lspk_IF=spk_IF{ii};  lspk_GLM=spk_GLM{ii};
            for tt=conds{ds,ntnum(ii)}
                spk=[lspk_act(tt,:) lspk_IF(tt,:) lspk_GLM(tt,:)];
                spk=cellfun(@(x) x(:)',spk,'uni',0);
                lenspk=cellfun(@length,spk);
                % exclude trial with too few spikes
                if(mean(cellfun(@length,spk))<spknumlim),continue;end
                % find ideal lag
                ldist=zeros(length(lags),14);
                for ll=1:length(lags)
                    s=[spk(1) cellfun(@(x) x+lags(ll),spk(2:end),'uni',0)];
                    dd=spkd_wrapper(s,2000);
                    ldist(ll,:)=dd(1,2:end);
                end
                [~,lagid]=min(ldist,[],1);
                spk(2:end)=cellfun(@(x,y) x+y,spk(2:end),num2cell(lags(lagid))','uni',0);
                % compute dist
                for qq=1:nq
                    distmatrix=spkd_wrapper(spk,q(qq))./bsxfun(@plus,lenspk,lenspk');
                    adIF=distmatrix(1:5,1:5);
                    ldActInt(tt,qq)=nanmean(column(adIF(~eye(5))));
                    ldIF(tt,qq)=nanmean(column(distmatrix(1:5,6:10)));
                    ldGLM(tt,qq)=nanmean(column(distmatrix(1:5,11:15)));
                end
            end
            dActInt{ii}=ldActInt;  dIF{ii}=ldIF; dGLM{ii}=ldGLM;
        end
        % save
        save(['mat/timing/spkd_' filename],'q','lags','spknumlim','distmatrix',...
            'dActInt','dIF','dGLM')
    else
        fprintf('reloading spike distances: '); tic
        load(['mat/timing/spkd_' filename])
    end
    toc
    
    %% same spike dist computation, but with jitter (only measured)
    if 1
        fprintf('computing spike distances with jitter: '); tic
        % select subset of conditions models have been fitted on
        conds={1:120 1:120 1:251; ... sines
            41:60 41:60 1:60; ... noises
            [1:55 76:80] [1:55 76:80] 1:138}; % diharms
        % parameters
        q=logspace(0.1,3.9,50); nq=length(q);
        lags=(-.005:.00025:.005)';
        spknumlim=3;
        %  jit=[0 .001 .002 .005 .01];
        jit=logspace(-4,-1.3,100)
        njit=length(jit);
        
        % compute
        dActIntJIT=cell(naff,njit);
        parfor jj=1:njit
            jj
            for ii=1:naff
                ldActInt=nan(ncond,nq);
                lspk_act=spk_act{ii};
                for tt=conds{ds,ntnum(ii)}
                    spk=lspk_act(tt,:);
                    spk=cellfun(@(x) x(:)',spk,'uni',0);
                    %add jitter
                    spk=[spk cellfun(@(x) x+randn(size(x))*jit(jj),spk,'uni',0)];
                    lenspk=cellfun(@length,spk);
                    % exclude trial with too few spikes
                    if(mean(cellfun(@length,spk))<spknumlim),continue;end
                    % find ideal lag
                    ldist=zeros(length(lags),9);
                    for ll=1:length(lags)
                        s=[spk(1) cellfun(@(x) x+lags(ll),spk(2:end),'uni',0)];
                        dd=spkd_wrapper(s,2000);
                        ldist(ll,:)=dd(1,2:end);
                    end
                    [~,lagid]=min(ldist,[],1);
                    spk(2:end)=cellfun(@(x,y) x+y,spk(2:end),num2cell(lags(lagid))','uni',0);
                    % compute dist
                    for qq=1:nq
                        distmatrix=spkd_wrapper(spk,q(qq))./bsxfun(@plus,lenspk,lenspk');
                        adIF=distmatrix(1:5,1:5);
                        ldActInt(tt,qq)=nanmean(column(distmatrix(1:5,6:10)));
                    end
                end
                dActIntJIT{ii,jj}=ldActInt;
            end
        end
    end
    toc
    
    %% plot ISI hist
    if 0
        % choose trial (61th is the sum of all trials)
        trial=size(isih_act,1);
        
        newfig('ISI','pos',[1 3 25 15])
        ax=zeros(naff,1);
        affidx=[0 0 0];
        goffset=[0 5 15];
        for ii=1:naff
            typ=ntnum(ii);
            affidx(typ)=affidx(typ)+1;
            ax(ii)=subplot(4,5,goffset(typ)+affidx(typ)); hold on;
            
            plot(isi_bin*1000,isih_act{trial,ii},'color',[1 .7 .7])
            hl(ii,1)=plot(isi_bin*1000,mean(isih_act{trial,ii},2),'r',...
                'displayname','Actual');
            
            plot(isi_bin*1000,isih_IF{trial,ii},'color',[.7 .7 1])
            hl(ii,2)=plot(isi_bin*1000,mean(isih_IF{trial,ii},2),'b',...
                'displayname','IF');
            
            plot(isi_bin*1000,isih_GLM{trial,ii},'color',[.4 .9 .4 .4])
            hl(ii,3)=plot(isi_bin*1000,mean(isih_GLM{trial,ii},2),'col',[0 .6 0],...
                'displayname','GLM');
            
            xlim(isi_bin([1 end])*1000)
        end
        
        set(ax,'xscale','log','xtick',[1 10 100 1000])
        xlabel(ax(14),'ISI [ms]')
        ylabel(ax(14),'Count')
        lll=lgd(hl(end,:)); set(lll,'pos',get(lll,'pos')+[.1 0 0 0])
        p=get(ax(1),'pos');
        axt(1)=axes('pos',[.03 p(2) .03 p(4)],'color',affcol(1));
        p1=get(ax(10),'pos'); p2=get(ax(5),'pos');
        axt(2)=axes('pos',[.03 p1(2) .03 p2(2)+p2(4)-p1(2)],'color',affcol(2));
        p=get(ax(14),'pos');
        axt(3)=axes('pos',[.03 p(2) .03 p(4)],'color',affcol(3));
        set(axt,'xcolor','none','ycolor','none');
        prop={'horiz','center','vert','middle','rot',90,'fontsi',16,'fontw','bold'};
        for ii=1:3
            textnorm(axt(ii),.5,.5,aff{ii},prop{:})
        end
        
        screenshot(['figs/timing/timing_' datasets{ds} '_isi_hist_n' num2str(trial)])
    end
    %% computing ISI histograms CORR for diff bin size
    if 0
        binl=[2 5 10 20 50 100 200 500];
        
        % corr value for rand var
        randcorr=zeros(length(binl),1);
        for jj=1:length(binl)
            adIF=corr(rand(binl(jj),100));
            randcorr(jj)=mean(abs(adIF(~eye(size(adIF)))));
        end
        
        newfig('ISI-hist-CORR-bin-size-effect','pos',[1 3 25 15])
        ax=zeros(naff,1);
        affidx=[0 0 0];
        goffset=[0 5 15];
        for ii=1:naff
            
            corr_IF=zeros(length(binl),ncond);
            corr_GLM=zeros(length(binl),ncond);
            
            typ=ntnum(ii);
            affidx(typ)=affidx(typ)+1;
            ax(ii)=subplot(4,5,goffset(typ)+affidx(typ)); hold on;
            for jj=1:length(binl)
                edges=logspace(-3.2,0,binl(jj));
                binwidth = diff(edges);
                edges = [-Inf, edges(1:end-1)+binwidth/2, Inf];
                
                lisi_act=cellfun(@diff,spk_act{ii},'uni',0);
                lisi_act=num2cell(cellfun(@(x) histcountsmex(x,edges),lisi_act,'uni',0),2);
                lisih_act=cellfun(@(x) cat(1,x{:})',lisi_act,'uni',0);
                lisih_act{end+1}=sum(cat(3,lisih_act{:}),3);
                
                lisi_IF=cellfun(@diff,spk_IF{ii},'uni',0);
                lisi_IF=num2cell(cellfun(@(x) histcountsmex(x,edges),lisi_IF,'uni',0),2);
                lisih_IF=cellfun(@(x) cat(1,x{:})',lisi_IF,'uni',0);
                lisih_IF{end+1}=sum(cat(3,lisih_IF{:}),3);
                
                lisi_GLM=cellfun(@diff,spk_GLM{ii},'uni',0);
                lisi_GLM=num2cell(cellfun(@(x) histcountsmex(x,edges),lisi_GLM,'uni',0),2);
                lisih_GLM=cellfun(@(x) cat(1,x{:})',lisi_GLM,'uni',0);
                lisih_GLM{end+1}=sum(cat(3,lisih_GLM{:}),3);
                
                for kk=1:ncond
                    xx=corr([lisih_act{kk} lisih_IF{kk} lisih_GLM{kk}]);
                    corr_IF(jj,kk)=mean(column(xx(1:5,6:10)));
                    corr_GLM(jj,kk)=mean(column(xx(1:5,11:15)));
                end
            end
            
            h(1,:)=plotstd(binl,corr_IF,@(x) deal(nanmean(x,2),nansem(x,1,2)),'col','b','displayname','IF');hold on
            h(2,:)=plotstd(binl,corr_GLM,@(x) deal(nanmean(x,2),nansem(x,1,2)),'col',[0 .6 0],'displayname','GLM');
            plot(binl,randcorr,'--','col',[.5 .5 .5]);
            set(gca,'xscale','log','xlim',binl([1 end]))
            grid
            
        end
        lll=lgd(h(1:2,:)); set(lll,'pos',get(lll,'pos')+[.1 0 0 0])
        
        xlabel(ax(14),'Number of bins')
        ylabel(ax(14),'Correlation')
        
        p=get(ax(1),'pos');
        axt(1)=axes('pos',[.03 p(2) .03 p(4)],'color',affcol(1));
        p1=get(ax(10),'pos'); p2=get(ax(5),'pos');
        axt(2)=axes('pos',[.03 p1(2) .03 p2(2)+p2(4)-p1(2)],'color',affcol(2));
        p=get(ax(14),'pos');
        axt(3)=axes('pos',[.03 p(2) .03 p(4)],'color',affcol(3));
        set(axt,'xcolor','none','ycolor','none');
        prop={'horiz','center','vert','middle','rot',90,'fontsi',16,'fontw','bold'};
        for ii=1:3
            textnorm(axt(ii),.5,.5,aff{ii},prop{:})
        end
        
        screenshot(['figs/timing/timing_' datasets{ds} '_isi_binning_size'])
    end
    
    %% spike dist
    if 0
        newfig('spike dist','pos',[1 3 25 15])
        ax=zeros(naff,1);
        affidx=[0 0 0];
        goffset=[0 5 15];
        yl=zeros(naff,2);
        for ii=1:naff
            typ=ntnum(ii);
            affidx(typ)=affidx(typ)+1;
            ax(ii)=subplot(4,5,goffset(typ)+affidx(typ)); hold on;
            clear h
            h(1,:)=plotci(1000./q,dIF{ii}',[],'col','b','displayn','IF'); hold on
            h(2,:)=plotci(1000./q,dGLM{ii}',[],'col',[0 .6 0],'displayn','LNP');
            h(3,:)=plotci(1000./q,-dIF{ii}'+dGLM{ii}',[],'col',[.8 .2 0],'displayn','diff: LNP - IF'); hold on;
            plot(1000./q,zeros(size(q)),'k--')
            set(gca,'xscale','log')
            grid; axis tight
            xlim(1000./q([end 1]))
            yl(ii,:)=ylim;
        end
        set(ax,'xtick',[1 10 100],'ytick',[0 1],'xticklabels',[1 10 100],...
            'ylim',[min(yl(:,1)) 1])
        lll=lgd(h(1:3,1)); set(lll,'pos',get(lll,'pos')+[.18 0 0 0])
        
        xlabel(ax(14),'Time scale [ms]')
        ylabel(ax(14),'Spike distance')
        
        p=get(ax(1),'pos');
        axt(1)=axes('pos',[.03 p(2) .03 p(4)],'color',affcol(1));
        p1=get(ax(10),'pos'); p2=get(ax(5),'pos');
        axt(2)=axes('pos',[.03 p1(2) .03 p2(2)+p2(4)-p1(2)],'color',affcol(2));
        p=get(ax(14),'pos');
        axt(3)=axes('pos',[.03 p(2) .03 p(4)],'color',affcol(3));
        set(axt,'xcolor','none','ycolor','none');
        prop={'horiz','center','vert','middle','rot',90,'fontsi',16,'fontw','bold'};
        for ii=1:3
            textnorm(axt(ii),.5,.5,aff{ii},prop{:})
        end
        screenshot(['figs/timing/timing_' datasets{ds} '_spike_dist'])
    end
    %% spike dist summary
    if 0
        newfig('spike dist 2','pos',[1 3 8.7 8])
        adIF=cellfun(@nanmedian,dIF,'uni',0);
        adIF=cat(1,adIF{:})';
        adGLM=cellfun(@nanmedian,dGLM,'uni',0);
        adGLM=cat(1,adGLM{:})';
        adiff=cellfun(@(x,y) nanmedian(x-y),dGLM,dIF,'uni',0);
        adiff=cat(1,adiff{:})';
        
        cols=[repmat([0 0 0],5,1);repmat([0 0 1],5,1);repmat([0 .6 0],5,1)];
        
        ax=subplot_ax(2,3,'nextplot','add'); clear yl
        nidx=[3 6 14];
        tr=[35 30 40];
        for ii=1:3
            spk=[spk_act{nidx(ii)}(tr(ii),:)...
                spk_IF{nidx(ii)}(tr(ii),:)...
                spk_GLM{nidx(ii)}(tr(ii),:)];
            spk=cellfun(@(x) x(:)',spk,'uni',0);
            ldist=zeros(length(lags),14);
            for ll=1:length(lags)
                s=[spk(1) cellfun(@(x) x+lags(ll),spk(2:end),'uni',0)];
                dd=spkd_wrapper(s,200);
                ldist(ll,:)=dd(1,2:end);
            end
            [~,lagid]=min(ldist,[],1);
            spk(2:end)=cellfun(@(x,y) x+y,spk(2:end),num2cell(lags(lagid))','uni',0);
            plot_spikes(spk,'color',cols,'parent',ax(ii))
            
            %set(ax(ii),'pos',get(ax(ii),'pos')+[-.02 .18 .00 -.23])
            id=ntnum==ii;
            
            plot(ax(ii+3),1000./q,zeros(size(q)),'-','col',[.4 .4 .4])
            plot(ax(ii+3),1000./q,adIF(:,id),'col',[0 0 1 .3])
            plot(ax(ii+3),1000./q,adGLM(:,id),'col',[0 .6 0 .3])
            plot(ax(ii+3),1000./q,adiff(:,id),'col',[.8 .2 0 .3])
            h(1)=plot(ax(ii+3),1000./q,mean(adIF(:,id),2),'col',[0 0 1],'displayname','IF');
            h(2)=plot(ax(ii+3),1000./q,mean(adGLM(:,id),2),'col',[0 .6 0],'displayname','LNP');
            h(3)=plot(ax(ii+3),1000./q,mean(adiff(:,id),2),'col',[.8 .2 0],'displayname','DIFF');
            axis(ax(ii+3),'tight'); set(ax(ii),'xlim',[.15 .25],'ylim',[-1 15])
            yl(ii,:)=ylim(ax(ii+3));
            title(ax(ii),aff{ii});
            if(ii==1)
                scalebar(ax(ii),'x','50 ms',.05,0,'y','',0,0)
            else
                scalebar(ax(ii),'x','',0,0,'y','',0,0)
            end
        end
        lll=lgd(h);set(lll,'pos',get(lll,'pos')+[.3 .12 0 0])
        xlabel(ax(5),'Time scale [ms]')
        
        ylabel(ax(4),'Norm spk dist')
        set(ax(4:6),'xscale','log','xlim',1000./q([end 1]),'ylim',[min(yl(:,1)) 1],...
            'xtick',[.99999 9.9999999 100],'xticklabels',[1 10 100],...
            'ytick',[0 1],'xgrid','on')
        
        screenshot(['figs/timing/timing_' datasets{ds} '_spike_dist_summary'],'png')
        screenshot(['figs/timing/timing_' datasets{ds} '_spike_dist_summary'],'pdf')
    end
    %% same for IF vs ActIntJIT
    if 0
        for jj=1:njit
            plotname=['spike_dist_act_IF_' datasets{ds} '_jit_' num2str(jit(jj)*1000) 'ms'];
            newfig(plotname,'pos',[1 3 25 15])
            ax=zeros(naff,1);
            affidx=[0 0 0];
            goffset=[0 5 15];
            yl=zeros(naff,2);
            
            dIFminusbaseline=cellfun(@(x,y) bsxfun(@plus,x,y(:,1)-x(:,1)),dIF,dActIntJIT(:,jj),'uni',0);
            dIFActsigntruc=cellfun(@(x,y) ttest(x,y,'Alpha',0.05,'tail','right'),dIFminusbaseline,dActIntJIT(:,jj),'uni',0);
            dIFActsign{jj}=cat(1,dIFActsigntruc{:});
            
            for ii=1:naff
                typ=ntnum(ii);
                affidx(typ)=affidx(typ)+1;
                ax(ii)=subplot(4,5,goffset(typ)+affidx(typ)); hold on;
                idx_sign_diff=find(dIFActsign{jj}(ii,:)>0,1);
                if isempty(idx_sign_diff), idx_sign_diff=length(q);end
                clear h
                h(1,:)=plotci(1000./q,dIF{ii}',[],'col','b','displayn','IF'); hold on
                h(2,:)=plotci(1000./q,dActIntJIT{ii,jj}',[],'col',[0 .6 0],'displayn','ActInt');
                %                 h(3,:)=plotci(1000./q,+dIF{ii}'-dActIntJIT{ii,jj}',[],'col',[.8 .2 0],'displayn','diff: IF - ActInt'); hold on;
                h(3,:)=plotstd(1000./q,+dIFminusbaseline{ii}'-dActIntJIT{ii,jj}',[],'col',[.8 .2 0],'displayn','diff: IF - ActInt'); hold on;
                plot(1000./q,zeros(size(q)),'k--')
                plot(1000./q([idx_sign_diff idx_sign_diff]),[0 1])
                set(gca,'xscale','log')
                grid; axis tight
                xlim(1000./q([end 1]))
                yl(ii,:)=ylim;
            end
            set(ax,'xtick',[1 10 100],'ytick',[0 .25 .5 .75 1],'xticklabels',[1 10 100],...
                'ylim',[mean(yl(:,1)) 1])
            lll=lgd(h(1:3,1)); set(lll,'pos',get(lll,'pos')+[.18 0 0 0])
            
            xlabel(ax(14),'Time scale [ms]')
            ylabel(ax(14),'Spike distance')
            
            p=get(ax(1),'pos');
            axt(1)=axes('pos',[.03 p(2) .03 p(4)],'color',affcol(1));
            p1=get(ax(10),'pos'); p2=get(ax(5),'pos');
            axt(2)=axes('pos',[.03 p1(2) .03 p2(2)+p2(4)-p1(2)],'color',affcol(2));
            p=get(ax(14),'pos');
            axt(3)=axes('pos',[.03 p(2) .03 p(4)],'color',affcol(3));
            set(axt,'xcolor','none','ycolor','none');
            prop={'horiz','center','vert','middle','rot',90,'fontsi',16,'fontw','bold'};
            for ii=1:3
                textnorm(axt(ii),.5,.5,aff{ii},prop{:})
            end
            screenshot(['figs/timing/timing_' plotname])
        end
        %%
        prctsign=cellfun(@(x) mean(x,2),dIFActsign,'uni',0);
        abovecrit=cat(2,prctsign{:})>=.5;
        result=[1000*jit;sum(abovecrit(ntnum==1,:));sum(abovecrit(ntnum==2,:));sum(abovecrit(ntnum==3,:))]';
        
    end
    %% spike dist summary WITH JITTER
    if 1
        only_sum_sum=1;
        cols=[repmat([0 0 0],5,1);repmat([0 0 1],5,1);repmat([0 .6 0],5,1)];
        if ~only_sum_sum
            fh1=newfig('spike dist w jitter summary','pos',[1 3 8.7 14]);
            ax=subplot_ax(4,3,'nextplot','add'); clear yl
        else
            fh2=newfig('spike dist w jitter summary summary','pos',[10 10 8.7 4.5]);
            ax2=subplot_ax(1,3,'nextplot','add');
        end
        
        if ~only_sum_sum
            [~,subjitidx]=min(abs(bsxfun(@minus,jit,[.001 .002 .005 .01]')),[],2);
            subjit=[.001 .002 .005 .01]; nsubjit=length(subjit);
            for jj=1:nsubjit
                adIF=cellfun(@nanmedian,dIF,'uni',0);
                adIF=cat(1,adIF{:})';
                adActIntJIT=cellfun(@nanmean,dActIntJIT(:,subjitidx(jj)),'uni',0);
                adActIntJIT=cat(1,adActIntJIT{:})';
                adiffmean=cellfun(@(x,y) nanmean(x-y),dIF,dActIntJIT(:,subjitidx(jj)),'uni',0);
                adiffmean=cat(1,adiffmean{:})';
                
                %adiff=cat(1,adiff{:})';
                
                for ii=1:3
                    id=ntnum==ii;
                    
                    plot(ax(ii+3*(jj-1)),1000./q,mean(adiffmean(1,id),2)*ones(size(q)),'-','col',[.95 .65 0])
                    plot(ax(ii+3*(jj-1)),1000./q(1./q>subjit(jj)),adIF(1./q>subjit(jj),id),'col',[0 0 1 .3])
                    plot(ax(ii+3*(jj-1)),1000./q(1./q>subjit(jj)),adActIntJIT(1./q>subjit(jj),id),'col',[0 .6 0 .3])
                    plot(ax(ii+3*(jj-1)),1000./q(1./q>subjit(jj)),adiffmean(1./q>subjit(jj),id),'col',[.8 .2 0 .3])
                    h(1)=plot(ax(ii+3*(jj-1)),1000./q(1./q>subjit(jj)),mean(adIF(1./q>subjit(jj),id),2),'col',[0 0 1],'displayname','simulated');
                    h(2)=plot(ax(ii+3*(jj-1)),1000./q(1./q>subjit(jj)),mean(adActIntJIT(1./q>subjit(jj),id),2),'col',[0 .6 0],'displayname','jittered');
                    h(3)=plot(ax(ii+3*(jj-1)),1000./q(1./q>subjit(jj)),mean(adiffmean(1./q>subjit(jj),id),2),'col',[.8 .2 0],'displayname','difference');
                    axis(ax(ii+3*(jj-1)),'tight');
                    yl(ii,:)=ylim(ax(ii+3*(jj-1)));
                    if jj==1
                        title(ax(ii),aff{ii},'pos',[10 1.25]);
                    end
                    
                    
                end
                text(10,1.1,['jitter=' num2str(subjit(jj)*1000) 'ms'],'horiz','center',...
                    'parent',ax(3*jj-1))
            end
        else
            for jj=1:njit
                adiff(:,jj)=cellfun(@(x,y) mean(bsxfun(@minus,x(:,1./q>jit(jj))-y(:,1./q>jit(jj)),x(:,1)-y(:,1)),2),dIF,dActIntJIT(:,jj),'uni',0);
            end
        end
        if ~only_sum_sum
            lll=lgd(h);set(lll,'pos',get(lll,'pos')+[.25 .04 0 0])
            xlabel(ax(11),'Time scale [ms]')
            ylabel(ax(10),'Normalized spike distance','horiz','left','pos',[0.01 min(yl(:,1))])
            set(ax,'xscale','log','xlim',1000./q([end 1]),'ylim',[min(yl(:,1)) 1],...
                'xtick',[.99999 9.9999999 100],'xticklabels',[1 10 100],...
                'ytick',[0 1],'YMinorTick','on','xgrid','on','ygrid','on')
            screenshot(fh1,['figs/timing/timing_' datasets{ds} '_spike_distwjit_summary'],'png')
            screenshot(fh1,['figs/timing/timing_' datasets{ds} '_spike_distwjit_summary'],'pdf')
        else
            for aa=1:3,
                affval{aa}=[];
                plot(ax2(aa),[1 1 nan 10 10],[-.097 .2 nan -.097 .2],'col',[.87 .87 .87],'linew',1.2)
                plot(ax2(aa),[.2 60],[0 0],'k')
            end
            for aa=1:naff
                aadiff{aa}=cat(2,adiff{aa,:});
                plot(ax2(ntnum(aa)),1000*jit',nanmean(aadiff{aa},1)','col',[.8 .2 0 .3])
                if ntnum(aa)==2,xlabel(ax2(ntnum(aa)),'Jitter SD [ms]'); end
                if ntnum(aa)==1,ylabel(ax2(ntnum(aa)),'Norm. dist. diff. [-]');end
                affval{ntnum(aa)}(:,end+1)=nanmean(aadiff{aa},1)';
                tresult(aa,:)=ttest(aadiff{aa},0,'tail','right');
                tval(aa)=jit(find(tresult(aa,:)==0,1)-1);
                plot(ax2(ntnum(aa)),1000*tval([aa aa]),[-.09 -.06],'k','linew',1)
            end
            for aa=1:3
                h=plot(ax2(aa),1000*jit',mean(affval{aa},2),'col',[.8 .2 0],...
                    'displayname','difference');
                title(ax2(aa),aff{aa})
            end
            set(ax2,'xtick',[1 10],'xticklabels',[1 10],'xlim',[.4 15],...
                'ygrid','on','ylim',[-.1 .2],'xscale','log')
            %lgd(h,'location','NW')
            text(70,0,'TouchSim precision is','horiz','center','vert','middle','rot',-90)
            text(30,0,'<- Worse   Better ->','horiz','center','vert','middle','rot',-90)
            screenshot(fh2,['figs/timing/timing_' datasets{ds} '_spike_distwjit_summary_diff'],'png')
            screenshot(fh2,['figs/timing/timing_' datasets{ds} '_spike_distwjit_summary_diff'],'pdf')
        end
    end
    
%     figure
%     imagesc(tresult)
    %%
    if 0
        newfig('test')
        for ii=1:50
            subplot(121)
            boxplot([dIFminusbaseline{9}(:,ii),dActInt{9}(:,ii)])
            subplot(122)
            plot(dIFminusbaseline{9}(:,ii),dActInt{9}(:,ii),'o')
            refline(1,0)
            nanmean([dIFminusbaseline{9}(:,ii),dActInt{9}(:,ii)])
            ttest2(dIFminusbaseline{9}(:,ii),dActInt{9}(:,ii),'Alpha',0.01)
            pause
        end
    end
    %% in the measured data spike dist
    if 0
        newfig('spike dist inter data','pos',[10 10 5 4])
        adActInt=cellfun(@nanmedian,dActInt,'uni',0);
        adActInt=cat(1,adActInt{:})';
        
        for ii=1:3
            id=ntnum==ii;
            
            plot(1000./q,adActInt(:,id),'col',[affcol(ii) .3]);hold on
            h=plot(1000./q,mean(adActInt(:,id),2),'col',affcol(ii),'displayname',aff{ii});
        end
        
        set(gca,'xscale','log',...
            'xlim',1000./q([end 1]),'ylim',[min(yl(:,1)) 1],...
            'xtick',[.99999 9.9999999 100],'xticklabels',[1 10 100],...
            'ytick',[0 1],'xgrid','on','box','off')
        xlabel('Time scale')
        ylabel('Norm. spike dist.')
        screenshot(['figs/timing/timing_' datasets{ds} '_spike_dist_actinter'],'pdf')
    end
    
    
    %% plot RASTERs
    if 0
        % choose freq
        %unique(stim.freq)
        f=300;
        
        trials=29:40;%find(stim.freq==f);
        [~,idx]=sort(stim.amp(trials),'descend');
        trials=trials(idx);
        
        newfig('RASTER','units','normalized','outerposition',[0 0 1 1])
        ax=zeros(naff,1);
        affidx=[0 0 0];
        goffset=[0 5 15];
        
        
        for ii=1:naff
            typ=ntnum(ii);
            affidx(typ)=affidx(typ)+1;
            ax(ii)=subplot(4,5,goffset(typ)+affidx(typ)); hold on;
            
            for jj=1:length(trials)
                %         if(length(cat(1,act_spikes{ii}{trials(jj),:}))<10)
                %             break
                %         end
                plot_spikes(spk_act{ii}(trials(jj),:),'parent',ax(ii),...
                    'color','r','hold','on','neuron_offset',0+(jj-1)*15);
                plot_spikes(spk_IF{ii}(trials(jj),:),'parent',ax(ii),...
                    'color','b','hold','on','neuron_offset',5+(jj-1)*15);
                plot_spikes(spk_GLM{ii}(trials(jj),:),'parent',ax(ii),...
                    'color',[0 .6 0],'hold','on','neuron_offset',10+(jj-1)*15);
            end
            
            xlim([0 stim.len(trials(1))])
            ylim([0 .001+(jj-1)*15])
        end
        
        axt(1)=axes('pos',get(ax(1),'pos')-[.07 0 .1 0],'color',affcol(1));
        p1=get(ax(10),'pos'); p2=get(ax(5),'pos');
        axt(2)=axes('pos',p1+[-.07 0 -.1 p2(2)+p2(4)-p1(2)-p1(4)],'color',affcol(2));
        axt(3)=axes('pos',get(ax(14),'pos')-[.07 0 .1 0],'color',affcol(3));
        set(axt,'xcolor','none','ycolor','none');
        aff={'SA1','RA','PC'};
        prop={'horiz','center','vert','middle','rot',90,'fontsi',16,'fontw','bold'};
        for ii=1:3
            textnorm(axt(ii),.5,.5,aff{ii},prop{:})
        end
        
        
        screenshot(['figs/timing/timing_' datasets{ds} '_raster'])
    end
    
    %% plot rasters
    if 0
        cols=lines(2);
        newfig(['rasters_' datasets{ds}],'pos',[1 3 25 15])
        ax=zeros(naff,1);
        for kk=1:138
            clf
            affidx=[0 0 0];
            goffset=[0 5 15];
            yl=zeros(naff,2);
            for ii=1:naff
                typ=ntnum(ii);
                affidx(typ)=affidx(typ)+1;
                ax(ii)=subplot(4,5,goffset(typ)+affidx(typ));
                plot_spikes([spk_act{ii}(kk,:) spk_IF{ii}(kk,:)],...
                    'col',cols([1 1 1 1 1 2 2 2 2 2],:),'par',ax(ii))
            end
            subplot(4,5,5)
            textnorm(gca,.5,.5,num2str(kk)); set(gca,'visible','off')
            pause
        end
    end
    
end 
    
    
