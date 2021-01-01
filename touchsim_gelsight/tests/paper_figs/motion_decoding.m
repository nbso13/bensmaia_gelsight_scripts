clear all;close all;clc

% in case you want to reload data 
% (if not, it will that long the compute responses)
flag_reload=0;

flag_video=0;
affcol=[affcol;.4 .4 .4];

if flag_reload
    % load ../mat/motion_decoding_16dirs_12_08_2016_12_16_52
    % load ../mat/motion_decoding_16dirs_11_08_2016_09_27_38
    % load ../mat/motion_decoding_8dirs_11_08_2016_10_00_02
    %load ../mat/motion_decoding_16dirs_15_08_2016_10_36_57
    load ../mat/motion_decoding/motion_decoding_16dirs_25_08_2016_11_13_14.mat
    %load ../mat/motion_decoding_16dirs_09_09_2016_13_03_51
    naff=length(ap.afferents);
    ndir=size(spk,2);
    nrep=size(spk,3);
else
    % EXPERIMENT configuration
    ndir=16;
    nreppertex=5;
    ntex=3;
    nrep=nreppertex*ntex;
    posjittersd=.5;
    rad=4;
    
    if(~flag_video)
        pos=-8:.01:8;
        sf=2000;
    else
        disp('warning: changed pos vector for video purpose')
        pos=-8:.1:8;
        sf=200;
    end
    
    % Afferent population
    fprintf('Affpop ');tic
    ap=affpop_hand('D2d',1,[],'delay',false);
    ap.afferents(hypot(ap.location(:,1),ap.location(:,2))>rad)=[];
    naff=length(ap.afferents);
    toc
end

if(~exist('ntex','var'))
    ntex=3;
    nreppertex=5;
end

% triangulation
types={'SA1','RA','PC','All'};
for ii=1:3
    x.(types{ii})=ap.location(ap.(['i' types{ii}]),1);
    y.(types{ii})=ap.location(ap.(['i' types{ii}]),2);
    tri=delaunay(x.(types{ii}),y.(types{ii}));
    sss=[tri(:,[1 2]);tri(:,[2 3]);tri(:,[1 3])];
    seg.(types{ii})=unique(sort(sss,2),'rows');
    nseg.(types{ii})=size(seg.(types{ii}),1);
    lseg.(types{ii})=hypot(diff(x.(types{ii})(seg.(types{ii})),1,2),...
        diff(y.(types{ii})(seg.(types{ii})),1,2));
end

% texture
load mat/resamp_profilometry
tex=tex(1:241,:,:);
tex=bsxfun(@plus,tex(:,:,[15 37 55])/1000,permute([1.2;.8;.8],[2 3 1]));
tex(:,:,1)=1.5*tex(:,:,1);
gridXY=-12:0.1:12;
[xtex,ytex]=meshgrid(gridXY);
for ii=1:ntex
    F{ii} = griddedInterpolant(xtex',ytex',tex(:,:,ii)');
end

% stimulus
[xs,ys]=meshgrid(-rad:.1:rad);
xs=xs(:); ys=ys(:); rs=hypot(xs,ys);
xs(rs>rad)=[];ys(rs>rad)=[];rs(rs>rad)=[];

fprintf('\nExp conditions: \n\n')
fprintf('NDIR: %i\n',ndir)
fprintf('NAFF: %i\n',naff)
fprintf('NREP: %i\n',nrep)
fprintf('POSJ: %i\n',posjittersd)
fprintf('PRAD: %i\n',rad)
fprintf('SFRQ: %i\n',sf)
fprintf('\n\n')

%% compute response
affpos=ap.location;
if ~flag_reload
    % 1. stimulus
    %ramp=stim_ramp(1,length(pos)/sf,[0 0],sf,.1);
    fprintf('Stim   ');tic
    for ii=1:length(F)
        % video purpose
        if(flag_video)
            v = VideoWriter(['../vids/stim_' num2str(ii) '.mp4'],'MPEG-4');
            v.Quality = 100;
            open(v);
            tri = delaunay(xs,ys);
            
            for dd=1:ndir
                alpha=dd/ndir*2*pi;
                h=trisurf(tri,xs*cos(alpha)-ys*sin(alpha),xs*sin(alpha)+ys*cos(alpha),[],nan(size(xs))); view(2); shading interp;
                set(gca,'visible','off');colormap gray;
                set(gcf,'color','w')
                
                trace=zeros(length(pos),length(xs));
                for jj=21:length(pos)-20
                    trace(jj,:)=sqrt(rad^2-rs.^2).*F{ii}(xs,ys+pos(jj));
                    set(h,'cdata',trace(jj,:));drawnow;
                    writeVideo(v,getframe(gcf))
                end
                
                %s(ii)=Stimulus(bsxfun(@times,trace,ramp.trace),[xs ys],sf);
                s(ii)=Stimulus(trace,[xs ys],sf);
            end
            close(v)
        else
            trace=zeros(length(pos),length(xs));
            for jj=1:length(pos)
                trace(jj,:)=sqrt(rad^2-rs.^2).*F{ii}(xs,ys+pos(jj));
            end
            %s(ii)=Stimulus(bsxfun(@times,trace,ramp.trace),[xs ys],sf);
            s(ii)=Stimulus(trace,[xs ys],sf);
            
        end
        
        
    end
    toc
    if(flag_video)
        return
    end
    
    % 2. response (rotate aff pop location to simulate stimulus rotation)
    posjitter=randn(ndir,2,nrep)*posjittersd;
    spk=cell(naff,ndir,nrep);
    for ii=1:ndir
        % aff pos
        theta=-ii/ndir*pi*2;
        rot=[cos(theta) sin(theta);-sin(theta) cos(theta)];
        locaffpos=affpos*rot;
        
        for rr=1:nrep
            ap.location=bsxfun(@plus,locaffpos,posjitter(ii,:,rr));
            % response
            fprintf('Resp   ');tic
            rc=ap.response(s(ceil(rr/nreppertex)));
            spk(:,ii,rr)={rc.responses(:).spikes};
            toc
        end
    end
    ap.location=affpos;
    
    % save
    datestr=char(datetime('now','TimeZone','local',...
        'Format','dd_MM_yyyy_HH_mm_ss'));
    save(['../mat/motion_decoding_' num2str(ndir) 'dirs_' datestr],...
        'spk','rad','posjitter','posjittersd','ap','pos','sf')
end

%% time jitter spikes ?
timejittersd=0; % in ms

timejitter=num2cell(repmat(randn(1,ndir,nrep)*timejittersd*.001,naff,1,1));
spk=cellfun(@plus,spk,timejitter,'uni',0);
spk=cellfun(@(x) x(x>0),spk,'uni',0);
rates=cellfun(@length,spk);

%% plot row spikes
if 1
    exspk=cell(3,4,3,3);
    dircol=[0 0 0; .5 .5 .5;0 0 0;.5 .5 .5];
    rates=cellfun(@length,spk);
    for jj=1:3
        if jj~=1
            fh=newfig(['raster_' types{jj}],'pos',[5+(jj-1)*5 12 3 6.25]);
        else
            fh=newfig(['raster_' types{jj}],'pos',[5+(jj-1)*5 12 4 6.25]);
        end
        for ii=1:4
            for kk=1:3
                idx=ap.(['i' types{jj}]);
                locspikes=cellfun(@(x) x-.12,spk(idx,ii*2,kk*5),'uni',0);
                [~,idx2]=sort(rates(idx,ii*2,1),'descend');
                plot_spikes(locspikes(idx2(9:11)),'color',dircol(ii,:),...
                    'neuron_o',12*(ii-1)+4*(kk-1),...
                    'linew',1,'hold','on');
                
                idx=9*(ii-1)+3*(kk-1)+1;
                exspk(:,ii,kk,jj)=locspikes(idx2(9:11));
            end
        end
        set(gca,'xlim',[-.005 .35],'ylim',[-1 47],'ycolor','none')
        prop={'fontsize',10};
        scalebar(gca,prop,'x','100ms',.1,0)
        if jj==1
            set(gca,'pos',[.22 .08 .76 .83])
            prop={'fontsize',8,'rot',90,'horiz','center'};
            text(-.036,1.5,'\_\_',prop{:})
            text(-.04,1.5,'T1',prop{:})
            text(-.036,5.5,'\_\_',prop{:})
            text(-.04,5.5,'T2',prop{:})
            text(-.036,9.5,'\_\_',prop{:})
            text(-.04,9.5,'T3',prop{:})
            
            text(-.081,5.5,'\_\_\_\_\_\_\_',prop{:})
            text(-.081,17.5,'\_\_\_\_\_\_\_',prop{:})
            text(-.081,29.5,'\_\_\_\_\_\_\_',prop{:})
            text(-.081,41.5,'\_\_\_\_\_\_\_',prop{:})
            
            prop={'fontsize',10,'rot',90,'horiz','center'};
            text(-.085,5.5,'45',prop{:})
            text(-.085,17.5,'90',prop{:})
            text(-.085,29.5,'135',prop{:})
            text(-.085,41.5,'180',prop{:})
        else
            set(gca,'pos',[.01 .08 .98 .83])
        end
        title(types{jj})
        screenshot(fh,['../figs/motion_dir_raster' types{jj}],'pdf')
    end
end



%% classification analysis (based on MFR and on spike distance)
twin=[.05 .1 .2 .4 .8]; % win duration
xtvec=[.05 .1 .2 .4 .8];
res=[0 500];

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
        spikesshort=cellfun(@(x) x(x<twin(tt)),spk(aff,:,:),'uni',0);
        
        % spike count classification
        
        spcount=cellfun(@length,spikesshort);
        spcount=reshape(spcount,size(spcount,1),[])';
        
        % LDA
        fprintf('spike count classification LDA');tic
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
        
        % spike distance classification LDA
        fprintf('spike distance classification LDA \n')
        for ii=1:nres
            fprintf('aff %i/%i win %i/%i res %i/%i ... ',aa,4,tt,ntwin,...
                ii,nres);tic;
            dist=zeros(ndir*nrep,ndir*nrep,length(aff));
            for kk=1:length(aff)
                dist(:,:,kk)=spkd_wrapper(spikesshort(kk,:),res(ii));
            end
            % recover 'spatial' nd-coordinates from distances
            dist=sqrt(sum(dist.^2,3));
            coord=distmat2coord(dist);
            allperfdist(tt,aa,ii,:)=LDAclassL1O(coord,dirs,reps);
            toc
            
            % NN
            fprintf('spike dist classification NN');tic
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
% newfig('LDA','pos',[380 250 750 750])
% ax=subplot_optimal(nres);
% for jj=1:nres
%     axes(ax(jj))
%     box off
%     hold on
%
%     for ii=1:3
%         h(ii)=plotstd(twin,squeeze(allperfdist(:,ii,jj,:)),[],...
%             'displayname',types{ii},'color',affcol(ii,:));
%         plotstd(twin,squeeze(allperfcount(:,ii,:)),[],...
%             'color',affcol(ii,:),'linestyle','--')
%     end
%     plot(twin,ones(size(twin))/ndir,'k--')
%     title(['Resolution: ' num2str(1000/res(jj),'%3.2f') ' ms'])
%     lgd(h,'location','NW')
%     xlabel('Window duration [ms]')
%     ylabel('Score')
%     set(gca,'xtick',xtvec,'xticklabels',xtvec*1000,'xscale','log',...
%         'xlim',[twin(1) twin(end)],'ylim',[0 1])
% end

%% plot classification scores
newfig('NN','pos',[380 250 750 750])
ax=subplot_optimal(nres);
for jj=1:nres
    axes(ax(jj))
    box off
    hold on
    
    for ii=1:4
        h(ii)=plotci(twin,squeeze(allperfdistNN(:,ii,jj,:)),[],...
            'displayname',types{ii},'color',affcol(ii,:));
        plotci(twin,squeeze(allperfcountNN(:,ii,:)),[],...
            'color',affcol(ii,:),'linestyle','--')
    end
    plot(twin,ones(size(twin))/ndir,'k--')
    title(['Resolution: ' num2str(1000/res(jj),'%3.2f') ' ms'])
    lgd(h,'location','NW')
    xlabel('Window duration [ms]')
    ylabel('Score')
    set(gca,'xtick',xtvec,'xticklabels',xtvec*1000,'xscale','log',...
        'xlim',[twin(1) twin(end)],'ylim',[0 .5])
end

%% compute corr
n=15;
bw=.007;


allperfncorr=zeros(length(twin),3,nrep);
allperfncorrNN=zeros(length(twin),3,nrep*(nrep-1));
for pp=1:4
    if(pp<4)
        aff=find(ap.(['i' types{pp}]));
    end
    for tt=1:length(twin)
        fprintf('Corr: aff %i/%i win %i/%i ',pp,4,tt,length(twin));tic
        if(pp<4)
            spkloc=cellfun(@(x) x(x<twin(tt)),spk(aff,:,:),'uni',0);
            ncorr{pp,tt}=zeros(2*n+1,nseg.(types{pp}),ndir,nrep); %
            
            for ii=1:nseg.(types{pp})
                for jj=1:ndir
                    for kk=1:nrep
                        loccorr=spike_corr(spkloc{seg.(types{pp})(ii,1),jj,kk},...
                            spkloc{seg.(types{pp})(ii,2),jj,kk},bw,n);
                        ncorr{pp,tt}(:,ii,jj,kk)=loccorr;
                    end
                end
            end
            toc
            
            % select the most spiking pairs
            [~,idx]=sort(nanmean(nanmean(nansum(ncorr{pp,tt},1),3),4),'descend');
            idx=idx(1:round(nseg.(types{pp})/8));
            %idx=idx(1:min(end,5000));
            ncorr{pp,tt}=reshape(ncorr{pp,tt}(:,idx,:,:),(2*n+1)*length(idx),[])'; % nseg.(types{pp})
            
            ncorrnext=ncorr{pp,tt};
        else
            ncorrnext=[ncorr{1,tt} ncorr{2,tt} ncorr{3,tt}];
            ncorrex=ncorr{2,tt};
            %[~,idx]=sort(mean(ncorrnext));
            %ncorrnext=ncorrnext(:,idx(1:end/3));
        end
        
        if(tt==length(twin))
            %figure('pos',[80+600*(pp-1) 80 560 800])
            %imagesc(ncorr)
        end
        
        % LDA
        fprintf('spike corr classification LDA');tic
        allperfncorr(tt,pp,:)=LDAclassL1O(ncorrnext,dirs,reps);
        toc
        % NN
        fprintf('spike corr classification NN');tic
        dist=distmat(ncorrnext);
        blocks=mat2cell(dist,ndir*ones(1,nrep),ndir*ones(1,nrep));
        select=diag(true(1,nrep));
        blocks=blocks(~select);
        [~,b]=min(cat(3,blocks{:}));
        score=sum(bsxfun(@eq,squeeze(b),(1:ndir)'))/ndir;
        allperfncorrNN(tt,pp,:)=score;
        toc
    end
end

%% PAPER
newfig('perf paper','pos',[15 12 8 6])
graycol=.5*ones(1,3);
box off
hold on
%grid on
clear h
for ii=1:4
    h(ii)=plotstd(twin,squeeze(allperfncorr(:,ii,:)),@(x) deal(nanmean(x,2),nansem(x,1,2)),...
        'displayname',types{ii},'color',affcol(ii,:));
end
xlim([twin(1) twin(end)])
ylim([0 1])
plot(twin,ones(size(twin))/ndir,'-','color',graycol,'linew',1)
text(twin(1),1/ndir,'chance ','color',graycol,'horiz','right')
lgd(h,'location','NW')
xlabel('Window duration [ms]')
ylabel('p(correct)')
set(gca,'xtick',xtvec,'xticklabels',xtvec*1000,'xscale','log')

screenshot figs/motion_decoding_perf pdf

%% plot corr class perf (MNF dashed)
newfig('perf Xcorr','defaultAxesColorOrder',affcol,...
    'pos',[35 12 7.5 9])
graycol=.6*ones(1,3);
box off
hold on
clear h
for ii=1:4
    h(ii)=plotstd(twin,squeeze(allperfcount(:,ii,:)),[],...
        'displayname',types{ii},'color',affcol(ii,:));
    plotstd(twin,squeeze(allperfdist(:,ii,1,:)),[],...
        'color',affcol(ii,:),'linest','--')
end
xlim([twin(1) twin(end)])
ylim([0 1])
plot(twin,ones(size(twin))/ndir,'-','color',graycol,'linew',1)
text(twin(1),1/ndir,'chance ','color',graycol,'horiz','right')
lgd(h,'location','N')
xlabel('Window duration [ms]')
ylabel('Score')
set(gca,'xtick',xtvec,'xticklabels',xtvec*1000,'xscale','log')





%% SUPPL
newfig('perf suppl','defaultAxesColorOrder',affcol,...
    'pos',[25 12 7.5 9])
graycol=.6*ones(1,3);
box off
hold on
clear h
for ii=1:4
    h(ii)=plot(twin,mean(squeeze(allperfdist(:,ii,1,:)),2),...
        'displayname',types{ii},'color',affcol(ii,:));
    plot(twin,mean(squeeze(allperfdist(:,ii,2,:)),2),...
        'color',affcol(ii,:),'linest',':')
end
xlim([twin(1) twin(end)])
ylim([0 1])
plot(twin,ones(size(twin))/ndir,'-','color',graycol,'linew',1)
text(twin(1),1/ndir,'chance ','color',graycol,'horiz','right')
lgd(h,'location','NW')
xlabel('Window duration [ms]')
ylabel('p(correct)')
set(gca,'xtick',xtvec,'xticklabels',xtvec*1000,'xscale','log')

screenshot figs/motion_decoding_perf_suppl pdf

%% figures for the paper
% finger and contact area
newfig('afferents_on_hand','pos',[35   12   2.5 3.3])

[o,t,ppm,regionprop]=plot_hand(nan);
r=[cos(-t) -sin(-t);sin(-t) cos(-t)];
plot_hand('region','D2d','scalebar',0,'fill',[.85 .85 .85]); hold on
ang=(pi/8:pi/8:2*pi)';
quiver(ppm*rad*.5*cos(ang+t)+o(1),ppm*rad*.5*sin(ang+t)+o(2),...
    ppm*1.25*rad*cos(ang+t),ppm*1.25*rad*sin(ang+t),0,'-','filled',...
    'linew',1.5,'maxheadsize',.3,'color','y')
ang=linspace(0,2*pi,200);
plot(ppm*(rad*cos(ang))+o(1),ppm*(rad*sin(ang))+o(2),'k','linew',1.5)
% h=flipud(get(gcf,'child'))
% for ii=1:3
%     l=get(h(ii),'child'); delete(l(1:end-2))
% end
screenshot figs/motion_decoding_ca pdf

%% afferents and pairs
newfig('afferents_on_hand','pos',[35   14   10 1.5])
ax=[axes('pos',[0 0.01 1/3 .98]) axes('pos',[1/3 0.01 1/3 .98])...
    axes('pos',[2/3 0.01 1/3 .98])];
set(ax,'nextplot','add')
ang=linspace(0,2*pi,200);
for ii=1:3
    s=seg.(types{ii});
    xloc=x.(types{ii})(s);
    yloc=y.(types{ii})(s);
    xya=[xloc(:,1) yloc(:,1)]*r;
    xyb=[xloc(:,2) yloc(:,2)]*r;
    xy1=xya+.3*(xyb-xya);
    xy2=xya+.7*(xyb-xya);
    plot(ax(ii),rad*cos(ang),rad*sin(ang),'k','linew',1); hold on
    plot(ax(ii),[xy1(:,1) xy2(:,1)]',[xy1(:,2) xy2(:,2)]',...
        [xya(:,1) xyb(:,1)]',[xya(:,2) xyb(:,2)]','.','col',affcol(ii,:),...
        'linew',.5,'markersize',5)
    set(ax(ii),'xcolor','none','ycolor','none')
    axis(ax(ii),'equal')
    %title(ax(ii),types{ii},'color',affcol(ii,:),'fontw','bold','fontsi',18)
end
screenshot figs/motion_decoding_ap pdf


%% stimulus
newfig('stimuli','pos',[10   10   7.5   3.3])
ax=subplot_ax(1,3,'nextplot','add');
ang=linspace(0,2*pi,200);
for ii=1:3
    set(ax(ii),'pos',[.1+.3*(ii-1) 0.025 .33 .95])
    imagesc(tex(:,40:200,ii),'par',ax(ii)); axis(ax(ii),'equal')
    plot(ax(ii),80+40*cos(ang),40+40*sin(ang),'w')
    plot(ax(ii),80+40*cos(ang),200+40*sin(ang),'w--')
    quiver(ax(ii),80,36,0,164,0,'linew',2,'color','y','MaxHeadSize',.4)
end
annotation('textbox',[.12,.85,.07,.15],'string','T1','margin',5,...
    'vert','middle','horiz','center','edgec','none','backg','w','fontsize',14)
annotation('textbox',[.42,.85,.07,.15],'string','T2','margin',5,...
    'vert','middle','horiz','center','edgec','none','backg','w','fontsize',14)
annotation('textbox',[.72,.85,.07,.15],'string','T3','margin',5,...
    'vert','middle','horiz','center','edgec','none','backg','w','fontsize',14)
colormap(gray)
set(ax,'visible','off')

screenshot figs/motion_decoding_stim pdf

% stim dirs
newfig('stimdir','pos',[1000   358   140   140])
dircol=[0 0 0;.6 0 0; .7 .7 .7; 1 0 0];
hold on
dirvec=(1:ndir)'/ndir*2*pi;

for ii=1:ndir
    if(ii==2||ii==4||ii==6||ii==8)
        annotation('arrow',.5+[.1 .4]*cos(dirvec(ii)),...
            .5+[.1 .4]*sin(dirvec(ii)),'color',dircol(ii/2,:),...
            'linew',3,'headw',10,'headl',10)
    else
        annotation('arrow',.5+[.1 .35]*cos(dirvec(ii)),...
            .5+[.1 .35]*sin(dirvec(ii)),'linew',1,'headw',5,'headl',5)
    end
end

annotation('textbox',[.88,.48,.05,.05],'string','0','edgecolor','none','vert','middle','horiz','center')
annotation('textbox',[.48,.93,.05,.05],'string','\pi/2','edgecolor','none','vert','middle','horiz','center')
annotation('textbox',[0,.48,.05,.05],'string','\pi','edgecolor','none','vert','middle','horiz','center')
annotation('textbox',[.48,.04,.05,.05],'string','3\pi/2','edgecolor','none','vert','middle','horiz','center')

axis equal
set(gca,'visible','off')
screenshot figs/motion_decoding_stimdir pdf

%%
fh=newfig('xcorr','pos',[15 12 8 2.8])
ax=subplot_ax(1,3,'nextplot','add');
phi=-90;
amp=25;

nn=1;

for ii=1:3
    imagesc(1000*(-n*bw:bw:n*bw),(1:ndir)/ndir*2*pi,...
        imfilter(ncorrex((1:16)+80*(ii-1),(nn-1)*31+(1:31)),fspecial('aver',2),'replicate'),'par',ax(ii))
    plot(ax(ii),amp*sin(((0:ndir+1)-11)/ndir*2*pi+phi/180*pi),(0:ndir+1)/ndir*2*pi,'col',[0 0 0 .3])
    set(ax(ii),'pos',get(ax(ii),'pos')+[.05 .30 .00 -.3])
end
xlabel(ax(2),'Time lag [ms]')
ylabel(ax(1),'Direction')
set(ax,'box','off','xtick',-50:50:50,'xlim',[-82 75],'ylim',[0 2*pi],...
    'ytick',[])
set(ax(1),'ytick',0:pi:2*pi,'yticklabel',{'0','\pi','2\pi'},'clim',[0 140])
set(ax(2:3),'ycolor','none')

annotation('textbox',[.166,.84,.06,.1],'string','T1','margin',0,...
    'vert','middle','horiz','center','edgec','none','backg','w')
annotation('textbox',[.444,.84,.06,.1],'string','T2','margin',0,...
    'vert','middle','horiz','center','edgec','none','backg','w')
annotation('textbox',[.722,.84,.06,.1],'string','T3','margin',0,...
    'vert','middle','horiz','center','edgec','none','backg','w')
% colormap jet
screenshot(fh,'figs/motion_decoding_xcorr','pdf')