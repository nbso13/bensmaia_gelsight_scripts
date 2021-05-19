clear all;close all;clc

% in case you want to reload data 
% (if not, it will that long the compute responses)
flag_reload=0;
affcol=[affcol;.4 .4 .4];

if flag_reload
    load ../mat/orientation_decoding_32dir_25_08_2016_13_27_47
    %load ../mat/orientation_decoding_16dir_25_08_2016_10_02_13
    %load ../mat/orientation_decoding_32dir_19_08_2016_18_36_36
    %load ../mat/orientation_decoding_64dir_21_08_2016_19_59_56
    
    naff=length(ap.afferents);
    ndir=size(spikes,2);
    nrep=size(spikes,3);
else
    % EXPERIMENT configuration
    ndir=32;
    nrep=10;
    posjittersd=0.5;
    
    bar_width=1.6; % mm
    rad=4; % mm
    
    depth=1; % mm
    
    % Afferent population
    fprintf('affpop...');tic
    ap=affpop_hand('D2d',1,[],'delay',false);
    naff=length(ap.afferents);
    toc
end

[xpin,ypin]=meshgrid(-rad:.1:rad); xpin=xpin(:); ypin=ypin(:);
r=hypot(xpin,ypin);
xpin(r>rad)=[]; ypin(r>rad)=[]; r(r>rad)=[];
types={'SA1','RA','PC','All'};
dircol=[0 0 0;.8 0 0];

fprintf('\nExp conditions: \n\n')
fprintf('NDIR: %i\n',ndir)
fprintf('NAFF: %i\n',naff)
fprintf('NREP: %i\n',nrep)
fprintf('POSJ: %i\n',posjittersd)
fprintf('PRAD: %i\n',rad)
fprintf('BARW: %i\n',bar_width)
fprintf('\n\n')

%% Stimulus and response
if ~flag_reload
    affpos=ap.location;
    spikes=cell(naff,ndir,nrep);
    ramp=stim_ramp(depth,.2);
    
    % stimulus
    fprintf('stimulus...');tic
    zpin=ones(size(xpin));
    idx=(round(abs(ypin)*10000)/10000<bar_width/2);
    zpin(~idx)=[];
    s=Stimulus(bsxfun(@times,ramp.trace,zpin'),[xpin(idx) ypin(idx)],...
        ramp.sampling_frequency);
    toc
    
    posjitter=randn(ndir,2,nrep)*posjittersd;
    for ii=1:ndir
        % aff pos
        theta=-ii/ndir*pi;
        rot=[cos(theta) sin(theta);-sin(theta) cos(theta)];
        locaffpos=affpos*rot;
        for rr=1:nrep
            ap.location=bsxfun(@plus,locaffpos,posjitter(ii,:,rr));
            % response
            fprintf('response: dir %i/%i, rep %i/%i  ',ii,ndir,rr,nrep);tic
            ap.location=locaffpos;
            rc=ap.response(s);
            spikes(:,ii,rr)={rc.responses(:).spikes};
            toc
        end
    end
    ap.location=affpos;
    
    % save
    datestr=char(datetime('now','TimeZone','local',...
        'Format','dd_MM_yyyy_HH_mm_ss'));
    save(['../mat/orientation_decoding_' num2str(ndir) 'dir_' datestr],...
        'spikes','ap','posjitter','posjittersd','bar_width','rad','ramp')
end



%% time jitter spikes ?
timejittersd=0; % in ms

timejitter=num2cell(repmat(randn(1,ndir,nrep)*timejittersd*.001,naff,1,1));
spikes=cellfun(@plus,spikes,timejitter,'uni',0);
spikes=cellfun(@(x) x(x>0),spikes,'uni',0);
rates=cellfun(@length,spikes);
firstspk=cellfun(@(x) x(1),spikes(~cellfun(@isempty,spikes)));
tfirstspk=min([min(firstspk(ap.iSA1)) min(firstspk(ap.iRA)) min(firstspk(ap.iPC))]);

ratedir=mean(rates,3);
mrates=mean(ratedir,2);
sf=ramp.sampling_frequency;



%% plot row spikes
if 1
    prop={'fontsize',10};
    for jj=1:3
        fh(jj)=newfig(['raster_' types{jj}],'pos',[5+(jj-1)*5 12 2.5 6.25]);
        get(fh(jj),'pos')
        idx=ap.(['i' types{jj}]);
        locspikes=spikes(idx,:,1);
        [~,idx2]=sort(mrates(idx),'descend');
        for ii=1:2
            plot_spikes(locspikes(:,ii*2),'color',dircol(ii,:),...
                'title',types{jj},'psth',.005,'neuron_o',1000,'par',fh(jj),...
                'linew',.8,'hold','on','psth_duration',.45)
            ax=get(fh(jj),'child'); 
            hold(ax(2),'on')
            
            plot_spikes(locspikes(idx2(1:30),ii*2),'color',dircol(ii,:),...
                'title',types{jj},'neuron_o',(ii-1)*32,'par',ax(1),...
                'linew',.8,'hold','on')
            
        end
        axh=get(gcf,'children');
        set(axh(1),'pos',get(axh(1),'pos')+[0.03 -.04 .33 .03])
        set(axh(2),'pos',get(axh(2),'pos')+[0.03 -.03 .33 .025])
        
        set(axh(1),'ylim',[-1 60],'xcolor','none','ycolor','none','xlim',[-.015 .45])
        %set(axh(1).Title,'color',affcol(jj,:),'fontsize',18)
        %scalebar(axh(1),prop,'x','50 ms',.05,0.025)
        
        axis(axh(2),'tight')
        yl=get(axh(2),'ylim');xl=get(axh(2),'xlim');
        set(axh(2),'ylim',[yl(1)-.1*(yl(2)-yl(1)) yl(2)])
        set(axh,'xlim',[-.015 .45])
        scalebar(axh(2),prop,'x','50 ms',.05,0.025,'y','500 spks/s',25,0.08)
        plot(axh(2),1/sf:1/sf:.4,ramp.trace*yl(2)*.5,'col',[.7 .7 .7])
        
        screenshot(fh(jj),['figs/bar_orientation_raster' types{jj}],'pdf')
    end  
end

%% another raster plot
if 0
    figure
    ax=subplot_optimal(ndir);
    for ii=1:ndir
        plot_spikes(spikes(:,ii,1),'par',ax(ii))
    end
end

%% classification analysis (based on MFR and on spike distance)
res=[0 500]; % cost to move a spike per unit time (seconds)
twin=[.001 .0025 .005 .01 .025 .05 .075 .35 .4 .5]; % win duration
xtvec=[.001 .002 .005 .01 .02 .05 .1 .2 .5];

nres=length(res);
ntwin=length(twin);

allperfcount=zeros(length(twin),4,nrep);
allperfdist=zeros(length(twin),4,nres,nrep);

allperfcountNN=zeros(length(twin),4,nrep*(nrep-1));
allperfdistNN=zeros(length(twin),4,nres,nrep*(nrep-1));

dirs=repmat(1:ndir,1,nrep)';
reps=column(repmat(1:nrep,ndir,1));

for aa=1:4
    if(aa<4)
        aff=find(ap.(['i' types{aa}]));
    else
        aff=1:length(ap.afferents);
    end
    
    for tt=1:ntwin
        spikesshort=cellfun(@(x) x(x<(twin(tt)+tfirstspk)),...
            spikes(aff,:,:),'uni',0);
        spcount=cellfun(@length,spikesshort);
        spcount=reshape(spcount,size(spcount,1),[])';
        
        % spike count classification
        fprintf('spike count classification ');tic
        allperfcount(tt,aa,:)=LDAclassL1O(spcount,dirs,reps);
        toc
        
        % NN
        fprintf('spike count classification NN');tic
        distall=distmat(spcount);
        blocks=mat2cell(distall,ndir*ones(1,nrep),ndir*ones(1,nrep));
        select=diag(true(1,nrep));
        blocks=blocks(~select);
        [~,b]=min(cat(3,blocks{:}));
        score=sum(bsxfun(@eq,squeeze(b),(1:ndir)'))/ndir;
        allperfcountNN(tt,aa,:)=score;
        toc
        
        % spike distance classification
        fprintf('spike distance classification \n')
        for ii=1:nres
            fprintf('aff %i/%i win %i/%i res %i/%i ... ',aa,4,tt,ntwin,...
                ii,nres);tic;
            dist=zeros(ndir*nrep,ndir*nrep,length(aff));
            for kk=1:length(aff)
                dist(:,:,kk)=spkd_wrapper(spikesshort(kk,:),res(ii));
            end
            % recover 'spatial' nd-coordinates from distances
            dist=sqrt(sum(dist.^2,3));
            if(~any(dist(:)))
                allperfdist(tt,aa,ii,:)=1/ndir;
            else
                coord=distmat2coord(dist);
                allperfdist(tt,aa,ii,:)=LDAclassL1O(coord,dirs,reps);
            end
            
            blocks=mat2cell(dist,ndir*ones(1,nrep),ndir*ones(1,nrep));
            select=diag(true(1,nrep));
            blocks=blocks(~select);
            [~,b]=min(cat(3,blocks{:}));
            score=sum(bsxfun(@eq,squeeze(b),(1:ndir)'))/ndir;
            allperfdistNN(tt,aa,ii,:)=score;
            toc
        end
    end
end

%% plot classification scores
% newfig('Classification performances LDA','pos',[380 250 750 750])
% ax=subplot_optimal(nres);
% graycol=.6*ones(1,3);
% for jj=1:nres
%     axes(ax(jj));    box off;    hold on
%     for ii=1:4
%         h(ii)=plotstd(twin,squeeze(allperfdist(:,ii,jj,:)),[],...
%             'displayname',types{ii},'color',affcol(ii,:));
%         plotstd(twin,squeeze(allperfcount(:,ii,:)),[],...
%             'color',affcol(ii,:),'linestyle','--')
%     end
%     plot([twin(1) twin(end)],ones(2,1)/ndir,'-',...
%         [.05 .05],[0 1],':',...
%         [.45 .45],[0 1],':','color',graycol,'linew',1)
%     lgd(h,'location','NW');
%     xlabel('Window duration [ms]')
%     ylabel('Score')
%     text(twin(1),1/ndir,'chance ','color',graycol,'horiz','right')
%     %title(['Resolution: ' num2str(1000/res(jj),'%3.2f') ' ms'])
%     set(gca,'xtick',xtvec,'xticklabels',xtvec*1000,'xscale','log',...
%         'xlim',[twin(1) twin(end)],'ylim',[0 1])
% end
% screenshot figs/orientation_decoding_perf_res pdf

%% plot classification scores with timing (suppl mat)
newfig('Classification performances NN','pos',[10 12 7.5 7.5])
%ax=subplot_optimal(nres-1);
graycol=.5*ones(1,3);
rectangle('position',[.05 0.005 .3 .995],'edgec','none','facec',[.92 .92 .92])
box off;    hold on
for ii=1:4
    h(ii)=plotstd(twin,squeeze(allperfdist(:,ii,2,:)),[],...
        'displayname',types{ii},'color',affcol(ii,:));
    plot(twin,mean(squeeze(allperfdist(:,ii,1,:)),2),':',...
        'color',affcol(ii,:))
end
plot([twin(1) twin(end)],ones(2,1)/ndir,'-','color',graycol,'linew',1)
lgd(h,'location','NW');
xlabel('Window duration [ms]')
ylabel('p(correct)')
text(.007,1.03,'ramp','horiz','cen','color',graycol)
text(.13,1.03,'hold','horiz','cen','color',graycol)
text(twin(1),1/ndir,'chance ','color',graycol,'horiz','right')
set(gca,'xtick',xtvec,'xticklabels',xtvec*1000,'xscale','log',...
    'xlim',[twin(1) twin(end)],'ylim',[0 1])

screenshot figs/orientation_decoding_perf_suppl pdf

%% perf fig for paper
newfig('Classification performances PAPER','pos',[15 12 8 6])
graycol=.5*ones(1,3);
box off;    hold on
rectangle('position',[.05 0.005 .3 .995],'edgec','none','facec',[.92 .92 .92])
for ii=1:4
    h(ii)=plotstd(twin,squeeze(allperfdist(:,ii,1,:)),...
        [],'displayname',types{ii},'color',affcol(ii,:),'linestyle','-');
    %     plot(twin,mean(squeeze(allperfdist(:,ii,3,:)),2),'--',...
    %         'color',affcol(ii,:));
end
lgd(h,'location','NW');

plot([twin(1) twin(end)],ones(2,1)/ndir,'--','color',graycol,'linew',1)
xlabel('Window duration [ms]')
h=ylabel('p(correct)');
text(twin(1),(1/ndir)*2,'chance ','color',graycol,'horiz','right','vert','mid')
text(.007,1.07,'ramp','horiz','cen','color',graycol)
text(.13,1.07,'hold','horiz','cen','color',graycol)
set(gca,'xtick',xtvec(1:3:end),'xticklabels',xtvec(1:3:end)*1000,'xscale','log',...
    'xlim',[twin(1) twin(end)],'ylim',[0 1],...
    'pos',get(gca,'pos')+[0 0 0.08 0])

screenshot figs/orientation_decoding_perf pdf

%% figures for the paper

fh(1)=newfig('afferents_on_hand','pos',[10   12   6   3]);
fh(2)=newfig('afferents_on_hand_resp','pos',[10   16   6   6]);
ax=[];
ax(1,:)=gridfig(1,3,'par',fh(1),'handletype','matrix','titleline',5);
ax(2:3,:)=gridfig(2,3,'par',fh(2),'handletype','matrix','titleline',2.5);

% afferents
for ii=1:3
    subap=AfferentPopulation;
    subap.afferents=ap.afferents(ap.(['i' types{ii}]));
    axes(ax(1,ii))
    plot(subap,[],'region','D2d','onehand',true,'scalebar',0);
    axes(ax(2,ii))
    plot(subap,[],'region','D2d','onehand',true,'rate',...
        ratedir(ap.(['i' types{ii}]),1)/...
        max(ratedir(ap.(['i' types{ii}]),1)),...
        'scalebar',0);
    axes(ax(3,ii))
    plot(subap,[],'region','D2d','onehand',true,'rate',...
        ratedir(ap.(['i' types{ii}]),15)/...
        max(ratedir(ap.(['i' types{ii}]),15)),...
        'scalebar',0);
end
screenshot(fh(1),'figs/orientation_decoding_apall','pdf')
screenshot(fh(2),'figs/orientation_decoding_apresp','pdf')

%% stimulus
fh=newfig('stimulus_on_hand0','pos',[6   10   1.7   6]);
ax=gridfig(2,1,'par',fh,'handletype','matrix','titleline',2.5);

axes(ax(1))
[ori,th,ppm]=plot_hand('region','D2d','scalebar',0);
rot=[cos(-th) -sin(-th);sin(-th) cos(-th)];
pin=[xpin ypin]*rot; pin=pin*ppm; pin=bsxfun(@plus,pin,ori);
idx=(round(abs(sind(0)*xpin+cosd(0)*ypin)*10000)/10000<bar_width/2);
plot(pin(idx,1),pin(idx,2),'.','col',dircol(2,:))

axes(ax(2))
[ori,th,ppm]=plot_hand('region','D2d','scalebar',0);
rot=[cos(-th) -sin(-th);sin(-th) cos(-th)];
pin=[xpin ypin]*rot; pin=pin*ppm; pin=bsxfun(@plus,pin,ori);
idx=(round(abs(sind(90)*xpin+cosd(90)*ypin)*10000)/10000<bar_width/2);
plot(pin(idx,1),pin(idx,2),'col',dircol(1,:))
screenshot figs/orientation_decoding_stimALL pdf

%% test
% distfromcount=distmat(spcount);
% distfromspike=sqrt(sum(dist.^2,3));
%
% plot(distfromcount(:),distfromspike(:),'.')