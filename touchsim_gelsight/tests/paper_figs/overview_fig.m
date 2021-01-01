% figures for figure 1: overview of the model
clear;clc;close all

flag_reload=1;

if ~flag_reload
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
    A = interp2(x,y,A,X,Y);
    A=1-A;
    
    %% Design stimulus
    zers=zeros(round(length(A)*3/4),length(A));
    smat=[zers;A';zers]*1.5+1.5;
    steps=size(smat,1)-size(smat,2);
    rad=5;
    [x,y]=meshgrid(-rad:1/pins_per_mm:rad);
    r=hypot(x,y);
    xy=[x(:)-5.5 y(:)+.7];
    xy(r(:)>rad,:)=[];
    shape=zeros(steps,size(xy,1));
    sf=500;
    dur=size(shape,1)/sf;
    w=length(x)/2;
    spanmm=1; % half span range of edge cutting
    for ii=1:steps
        local=smat(ii:ii+length(A)-1,:);
        local=local.*sqrt(rad.^2-r.^2)/5;
        local(r>rad)=[];
        shape(ii,:)=local(:)';
    end
    
    shape=bsxfun(@times,shape,1+.1*sin(2*pi*40*(1:steps)/steps)');
    
    %% Create stim
    disp('Create stim')
    tic
    s=Stimulus(shape,xy,sf);
    toc
    
    %% Aff pop
    disp('Create Aff pop')
    tic
    ap=affpop_hand('D2');
    toc
    %%
    disp('Comp resp')
    tic
    r=ap.response(s);
    toc
    
    save ../mat/overview_resp s ap r
else
    load ../mat/overview_resp
end



fr=60;

%% PROFILES

figure(2);set(2,'unit','centimeter','pos',[5   10   3.2   4])
plot(s,'spatial','profile',fr); xlim([100 160]);ylim([420 488])
set(gca,'pos',[-.17 0 1.25 1])
delete(findobj(gcf,'type','UIControl'));
t=get(gca,'children');
delete(t(end))
axis off
figure(3);set(3,'unit','centimeter','pos',[10   10   3.2   4])
plot(s,'spatial','profiledyn',fr); xlim([100 160]);ylim([420 488])
set(gca,'pos',[-.17 0 1.25 1])
delete(findobj(gcf,'type','UIControl'));
t=get(gca,'children');
delete(t(end))
axis off

screenshot(2,'figs/methods_profile','png')
screenshot(3,'figs/methods_profiledyn','png')


%% STIMULUS
close all
a1=find(s.location(:,1)==-5.5-3&s.location(:,2)==.7);
a2=find(s.location(:,1)==-5.5&s.location(:,2)==.7);
a3=find(s.location(:,1)==-5.5+3&s.location(:,2)==.7);
idx=[a1 a2 a3];

figure(4); set(4,'unit','centimeter','pos',[5   10   3 3])
scatter(s.location(:,2),-s.location(:,1),1,s.trace(60,:),'filled');
axis equal; hold on
plot(s.location(idx,2),-s.location(idx,1),'ro','markersize',5)
xlim([-4.3 5.7])
ylim([.5 10.5])
set(gca,'visible','off','pos',[0.1 0.05 .8 1])
h=scalebar('x','2mm',2,0);
set(h,'linew',1.5)
set(h(2),'fontsize',10)

figure(5); set(5,'unit','centimeter','pos',[10   10   3 3])
mt=size(s.trace,1);
for ii=1:3
    ax(ii)=subplot(3,1,ii);
    set(ax(ii),'pos',get(ax(ii),'pos')+[-.04 .08 .15 .08])
    plot((1:mt)/mt,s.trace(:,idx(ii)),'r',[60 60]/mt,[0 4],'k:','linew',1);
    box off
    if ii==3
        h=scalebar('x','200ms',.2,.01,'y','2mm',2,.3);
        set(h,'linew',1.5)
        set(h([2 4]),'fontsize',10)
    end
    set(gca,'xlim',[-.025 1],'ylim',[0 4],'xcolor','none','ycolor','none')
end

screenshot(4,'figs/methods_stim1','png')
screenshot(5,'figs/methods_stim2','pdf')

%% collect strains and dyn part

allresp=[r.responses(:).propagated_struct];
strain=[allresp(:).stat_comp];
udyn=[allresp(:).dyn_comp];

figure(6);set(6,'pos',[400   100   160   530])
scatter(ap.location(ap.iRA,2),-ap.location(ap.iRA,1),[],strain(60,ap.iRA),'filled');
set(gca,'pos',[-.4 -.05 1.8 1.1])
axis off
figure(7);set(7,'pos',[800   100   160   530])
scatter(ap.location(ap.iRA,2),-ap.location(ap.iRA,1),[],log(udyn(60,ap.iRA)),'filled');
set(gca,'pos',[-.4 -.05 1.8 1.1])
axis off
screenshot(6,'figs/methods_prop1','pdf')
screenshot(7,'figs/methods_prop2','pdf')

%% affpop
figure(20)
set(20,'pos',[520 380 300 360])
plot(ap,[],'region','D2','scalebar',0,'rotate',-10)
ax=get(20,'children');
set(ax,'color','none')
axis(ax,'tight')
set(ax(3),'pos',[0+.02 .07 1/3 .91])
set(ax(2),'pos',[1/3 .07 1/3 .91])
set(ax(1),'pos',[2/3-.02 .07 1/3 .91])
a(1)=annotation('textbox','pos',[0 0 1/3 .1],'string','SA1');
a(2)=annotation('textbox','pos',[1/3 0 1/3 .1],'string','RA');
a(3)=annotation('textbox','pos',[2/3 0 1/3 .1],'string','PC');
set(a,'edgecolor','none','fonts',20)
screenshot(20,'figs/affpop','pdf')

%% sort and plot rasters
close all
notresp=r.rate<10;
r.affpop.afferents(notresp)=[];
r.responses(notresp)=[];
%r.responses=fliplr(r.responses);
%r.affpop.afferents=fliplr(r.affpop.afferents);

a=find(r.affpop.iSA1);
b=find(r.affpop.iRA);
c=find(r.affpop.iPC);

figure(8);set(8,'pos',[690   50   70   240])
plot(r,'linewidth',.5); hold on
set(gca,'box','off','xcolor','none','ycolor','none','xlim',[-.02 .311],...
    'ylim',[0 max(c)])
set(gca,'pos',get(gca,'pos')+[-.46 -.18 .39 .24])
text(-.11,mean(c),'PC','color',affcol(3),'rot',90,'vert','top','hor','center')
text(-.11,mean(b),'RA','color',affcol(2),'rot',90,'vert','top','hor','center')
text(-.11,mean(a),'SA1','color',affcol(1),'rot',90,'vert','top','hor','center')
screenshot(8,'figs/methods_raster','pdf')

%% p,v,a
close all
x=linspace(-1.5,1.5,200)';
y=-sin(2*x)+.75*sin(5*x);
y(:,2)=derive(y(:,1),1);
y(:,3)=derive(y(:,2),1);
yp=y; yp(yp<0)=0;
ym=y; ym(ym>0)=0; ym=-ym;
y=[y yp ym];
names={'p','v','a','pp','vp','ap','pm','vm','am'};

for ii=1:size(y,2)
    figure(8+ii); set(8+ii,'pos',[680   750   150   120])
    plot(x,y(:,ii),'col',[.6 .6 .6],'linew',4);
    set(gca,'ycolor','none','xtick',[],'box','off','linew',4,...
        'XAxisLocation','origin','pos',[0.05 0.05 .9 .9])
    screenshot(8+ii,['figs/methods_wave' names{ii}],'pdf')
end

%% lowpass
figure(18)
x=linspace(-5,5,200)';
y=zeros(size(x));
y(1:end/2)=1;
y(end/2+1:end)=1-x(end/2+1:end);
plot(x,y,'col',[.6 .6 .6],'linew',4); hold on
plot([0 0],[-4 3],'k--','linew',4)
set(18,'pos',[680   750   150   120])

set(gca,'linew',4,'xtick',[],'ytick',[],'box','off','ylim',[-4 3])
screenshot(18,'figs/methods_lowpass','pdf')

%% saturation
figure(19)
p=2;
x=linspace(-5,5,200)';
y=p*x./(p+abs(x));
plot(x,y,'col',[.6 .6 .6],'linew',4)
set(19,'pos',[680   750   150   120])
set(gca,'linew',4,'xtick',[],'ytick',[],'box','off',...
    'XAxisLocation','origin','YAxisLocation','origin')
screenshot(19,'figs/methods_sat','pdf')

