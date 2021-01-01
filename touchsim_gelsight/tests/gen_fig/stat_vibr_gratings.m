function stat_vibr_gratings()
%clear all; close all; clc

f=5000;
freqs=[0 5 10 20 40 80];
amps=[.325 .232 .169 .124 .045 .047];


% different frequencies
stat=0.325/700:0.325/700:0.325;
vib{1}=[stat(:);amps(1)*ones(19860,1);flipud(stat(:))];
for ff=2:6
    w=2*pi*freqs(ff);   t=1/f:1/f:10/freqs(ff);         vib{ff}=sin(w*t(:))*amps(ff)+1;
end

% general 3 periods gratings
gperiod=[1 2 4 6 8 10];
pinpermm=10;
for gg=1:6
    n=5*gperiod(gg);
    img{gg} = repmat([nan(1,n) 0*ones(1,n) nan(1,n) 0*ones(1,n) nan(1,n)]',1,3);
    shape{gg} = img2shape(img{gg},[0 NaN],pinpermm);
end
l=size(shape,1);

tic
s(length(freqs),length(gperiod))=Stimulus();

for ff=1:size(s,1)
    for gg=1:size(s,2)
        trace=[vib{ff}];
        t=(1:1:length(trace))/f;
        s(ff,gg) = stim_indent_shape(shape{gg},trace,[],f);
        fprintf('.')
    end
end
toc
fprintf('Stimulus processed.\n')


%%
model=[1:4 1:9];
class=[repmat({'SA1'},1,4) repmat({'RA'},1,9)];
nclass=length(class);
for ss=1:length(s(:))
    tic
    fprintf('S%02d/36',ss);
    ap=AfferentPopulation;
    for cc=1:nclass
        aploc=affpop_linear(0.5,max(s(ss).location(:,1)),class{cc},model(cc));
        ap.afferents=[ap.afferents aploc.afferents];
    end
    naff=size(ap.location,1);
    rc(ss)=ap.response(s(ss));
    spiketimes{ss}={rc(ss).responses.spikes};
    spikecount{ss}=cellfun(@length,{rc(ss).responses.spikes});
    srate{ss}=spikecount{ss}/s(ss).duration;
    fprintf('.')
    afflocs{ss}=ap.location;
    toc
end

fprintf('Response processed.\n')

%%

close all
for ii=1:36
    SArate{ii}=srate{ii}(1:end*4/13);
    RArate{ii}=srate{ii}(end*4/13+1:end);
    RAmeanr{ii}=nanmean(reshape(RArate{ii},[],9)',1);
    SAmeanr{ii}=nanmean(reshape(SArate{ii},[],4)',1);
end

fh=figure('position',[50 50 1500 600]);
set(fh,'defaultlinelinewidth',2,'DefaultAxesFontSize',14,'color','w')
col=get(0,'DefaultAxesColorOrder');
ax=gridfig(1,2,'parent',fh,'spacingh',2,...
    'marginb',0.1,'titleline',1);

yallshapes=[];
for ii=1:6
    shapeplot{ii}=img{ii}(end/5+1:3*end/5,1);
    yallshapes=[yallshapes;shapeplot{ii}];
end
yallshapes(yallshapes==0)=1;
yallshapes(isnan(yallshapes))=0;
xallshapes=((1:length(yallshapes))-1)/pinpermm;
initial_offset=cumsum([0 cellfun(@length,shapeplot)/pinpermm]);

for gg=1:6
    plot(ax(1),xallshapes,yallshapes,'col',[.6 .6 .6])
    plot(ax(2),xallshapes,yallshapes,'col',[.6 .6 .6])
    for ff=1:6
        idx=(gg-1)*6+ff;
        xaff=unique(afflocs{idx}(:,1));
        plot(ax(2),xaff(1:2*end/3)+initial_offset(gg),RAmeanr{idx}(1:2*end/3)/mean(RAmeanr{idx}),'.-','color',col(ff,:))
        plot(ax(1),xaff(1:2*end/3)+initial_offset(gg),SAmeanr{idx}(1:2*end/3)/mean(SAmeanr{idx}),'.-','color',col(ff,:))
    end
    %set(ax,'xlim',[0 x(end)])
end

legend(ax(2),'Grating','Static','5Hz','10Hz','20Hz','40Hz','80Hz')
set(ax,'ylim',[0 2.3])

xlabel(ax(1),'Position [mm]')
xlabel(ax(2),'Position [mm]')
ylabel(ax(1),'Norm. firing rate')
title(ax(1),'SA1','fontsize',18)
title(ax(2),'RA','fontsize',18)

% export_fig('stat_vibr_gratings01_curr.eps','-eps','-p0.06',gcf)
% export_fig('stat_vibr_gratings01_curr.png','-png',gcf)
screenshot(fh,'figs/stat_vibr_gratings01_curr','pdf')

%%

% RA=reshape(cellfun(@mean,RAmeanr'),6,6);
% SA=reshape(cellfun(@mean,SAmeanr'),6,6);
% RAv=reshape(cellfun(@std,RAmeanr'),6,6);
% SAv=reshape(cellfun(@std,SAmeanr'),6,6);
% xc=bsxfun(@plus,(1:6)',-.25:.1:.25);
% ax=gridfig(2,1,'spacing',5);
% errorbar(ax(1),xc,SA,SAv,'.-');hold on;
% errorbar(ax(2),xc,RA,RAv,'.-');hold on;
% set(ax,'xtick',1:6,'xticklabel',{'Static',5,10,20,40,80})
% ylabel(ax(1),'SA')
% ylabel(ax(2),'RA')
