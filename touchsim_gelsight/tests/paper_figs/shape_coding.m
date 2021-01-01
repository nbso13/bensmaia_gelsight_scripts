% shape_coding.m

try
    load shape_coding_stim.mat
    fprintf('Stimulus loaded.\n')
catch
    
    % Produce the weird dot pattern (Blake)
    smat=zeros(221,2100,'single');
    central=12:50:2100;
    interm=37:50:2100;
    
    smat(111,central)=1;
    smat(111+50,central)=1;
    smat(111-50,central)=1;
    smat(111+100,central)=1;
    smat(111-100,central)=1;
    
    smat(111+25,interm)=1;
    smat(111-25,interm)=1;
    smat(111+75,interm)=1;
    smat(111-75,interm)=1;
    
    diams=linspace(2.5,25,2100/25);
    
    dotmat=zeros(25,'single');
    [x,y]=meshgrid(-12:12,-12:12);
    r=hypot(x,y);
    for ii=1:2100/25
        dotmat=zeros(25,'single');
        dotmat(r<=diams(ii)/2)=.25;
        dotmat(r<=diams(ii)/2-.2)=.5;
        smat(:,(ii-1)*25+1:ii*25)=conv2(smat(:,(ii-1)*25+1:ii*25),dotmat,'same');
    end
    
    % compute shapes
    rad=5;
    [x,y]=meshgrid(-rad:.1:rad);
    r=hypot(x,y);
    xy=[x(:) y(:)];
    
    shape=zeros(round((2100-length(x))),size(xy,1));
    sf=500;
    dur=size(shape,1)/sf;
    w=(length(x)-1)/2;
    spanmm=1; % half span range of edge cutting
    for ii=1:size(shape,1)
        local=smat(111-w:111+w,(ii-1)+1:(ii-1)+length(x));
        shape(ii,:)=local(:)';
    end
    
    % Compute stimulus
    sf=[400,400,400];
    indmult=[280,370,620]/500;
    
    fprintf('Processing stimulus...\n')
    tic
    for ii=1:length(sf);
        fprintf('.');
        s(ii)=Stimulus(shape*indmult(ii),xy,sf(ii));
    end
    fprintf('Done   '),toc
    save shape_coding_stim.mat s w y sf smat shape
end

%% Compute response

ap=AfferentPopulation;

SAv=1:4;
RAv=1:9;
PCv=1:4;

naff=60;
for ii=0:naff, add_afferents(ap,'SA1',[0 ii/5-naff/10],'idx',SAv(ceil(rand*length(SAv))));end
for ii=0:naff, add_afferents(ap,'RA',[0 ii/5-naff/10],'idx',RAv(ceil(rand*length(RAv))));end
for ii=0:naff, add_afferents(ap,'PC',[0 ii/5-naff/10],'idx',PCv(ceil(rand*length(PCv))));end


disp('Processing Response')
tic;
for jj=1:3
    fprintf('.');
    rc(jj)=ap.response(s(jj));
end
fprintf('Done   '),toc

% Sort spikes
nrec=length(ap.afferents)/3;
spikes_t=cell(nrec,3,3);
for ii=1:length(rc)
    spikes_t(:,1,ii)={rc(ii).responses(ap.iSA1).spikes};
    spikes_t(:,2,ii)={rc(ii).responses(ap.iRA).spikes};
    spikes_t(:,3,ii)={rc(ii).responses(ap.iPC).spikes};
end

%% make figures

fh=figure('units','normalized','outerposition',[0 0 1 1],'DefaultAxesFontSize',14);
ax=gridfig(9,1,'parent',fh,'marginl',0.08,'marginb',.12,'spacingv',.2,'titleline',.1);
pstim=-w*.1:.1:size(shape,1)*.1+w*.1;
pcolor(ax(1),pstim(51:end-50),y(:,1),double(smat(111-w:111+w,51:end-50)));

delete(ax([2 6]))
ax=ax([1 3:5 7:9]);

for ii=1:3
    for jj=1:2
        plot_spikes(spikes_t(:,jj,ii),'par',ax(ii+3*(jj-1)+1),...
            'time_s',sf(ii)/10,'neuron_s',1/5,'neuron_o',-30.5);
    end
end

for ii=1:7
    shading(ax(ii),'flat')
end

linkaxes(ax(:),'xy')
%labels={'Stim.','280\mum','370um','620um','280um','370um','620um'};

for ii=1:7
    %axes(ax(ii)),text(202,0,labels{ii},'fontsize',24)
    xlabel(ax(ii),'');
    ylabel(ax(ii),'');
    set(ax(ii),'pos',get(ax(ii),'pos')+[-.05 0 0 0])
end

xlabel(ax(end),'Distance [mm]','FontSize',24)
ylabel(ax(3),'Distance [mm]','FontSize',24)
ylabel(ax(6),'Distance [mm]','FontSize',24)

set(fh,'color','w')
set(ax,'linewidth',1,'ylim',[-5 5],'tickdir','out',...
    'ticklength',[.005 .005],'xlim',[0 200],'box','on','fontsize',24)
title(ax(2),'SA1','fontsize',38)
title(ax(5),'RA','fontsize',38)
set(ax([1 2 4 5 7]),'yticklabel',{})
set(ax(2:end),'clim',[-100 1200])
set(ax([2 3 5 6]),'xtick',[])
set(ax([2 3 5 6]),'xticklabel',[])
set(ax(1:7),'xticklabel',[])
hotmap=hot(300);
colormap(fh,flipud(hotmap(end/5:end,:)))

%print(fh,sprintf('figs/phillips_blake%02d_curr.png',protocol),'-dpng')
screenshot(fh,'shape_coding1')
close(fh)
