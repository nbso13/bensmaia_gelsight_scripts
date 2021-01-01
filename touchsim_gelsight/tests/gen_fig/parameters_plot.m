function parameters_plot()
clear
close all
%clc

parameters=IF_parameters(1);

sap=reshape([parameters.sa{:}],length(parameters.sa{1}),[])';
rap=reshape([parameters.ra{:}],length(parameters.ra{1}),[])';
pcp=reshape([parameters.pc{:}],length(parameters.pc{1}),[])';

nsa = length(parameters.sa);
nra = length(parameters.ra);
npc = length(parameters.pc);



fig=figure('units','centimeters','pos',[25 9 17.8 9]);
ax=[subplot(2,11,1:2) subplot(2,11,3:5) subplot(2,11,6:8) subplot(2,11,9:11)...
    subplot(2,11,12:13) subplot(2,11,14:15) subplot(2,11,16:17) ...
    subplot(2,11,18:20) subplot(2,11,21:22)];

for ii=1:length(ax)
    set(ax(ii),'pos',get(ax(ii),'pos').*[1.25 1.05 .7 .8]+[-.08 .03 0 0])
end
set(ax(2),'pos',get(ax(2),'pos')+[0.025 0 0 0])
set(ax(3),'pos',get(ax(3),'pos')+[0.0125 0 0 0])
set(ax(8),'pos',get(ax(8),'pos')+[0.035 0 0 0])

set(ax,'nextplot','add','colororder',affcol)
dva=[2 4 6];
titles={'(2-3)','(4-5)','(6-7)'};
labels={'Quasi-static','Dynamic','Dynamic (d/dt)'};

%%
% if 0
%     % pos vel acc coef plot
%     
%     h=[];
%     for ii=1:3
%         h(end+1)=plot(ax(ii+1),(sap(:,dva(ii))),(sap(:,dva(ii)+1)),'.','color',affcol(1),'displayname','SA1');
%         h(end+1)=plot(ax(ii+1),(rap(:,dva(ii))),(rap(:,dva(ii)+1)),'.','color',affcol(2),'displayname','RA');
%         h(end+1)=plot(ax(ii+1),(pcp(:,dva(ii))),(pcp(:,dva(ii)+1)),'.','color',affcol(3),'displayname','PC');
%         
%         if(ii==2)
%             axx=axes('pos',[.65 .80 .06 .12])
%             h(end+1)=plot(sap(:,dva(ii)),sap(:,dva(ii)+1),'.','color',affcol(1)); hold on
%             h(end+1)=plot(rap(:,dva(ii)),rap(:,dva(ii)+1),'.','color',affcol(2)); box off
%             %set(axx,'xlim',[-.25 .25],'ylim',[-.25 .25])
%         end
%         
%         if(ii==3)
%             axx=axes('pos',[.88 .80 .06 .12])
%             h(end+1)=plot(rap(:,dva(ii)),rap(:,dva(ii)+1),'.','color',affcol(2)); box off
%             %set(axx,'xlim',[-.25 .25],'ylim',[-.25 .25])
%         end
%         
%         xlabel(ax(ii+1),[labels{ii} ' (+)'])
%         ylabel(ax(ii+1),[labels{ii} ' (-)'])
%         title(ax(ii+1),titles{ii})
%     end
%     %set(ax(2:4),'xscale','log','yscale','log')
%     set(ax(2),'xlim',[-.5 1.5],'ylim',[-.5 1.5])
%     set(ax(3),'xlim',[-280 280],'ylim',[-280 280])
%     set(ax(4),'xlim',[-250 1500],'ylim',[-250 1500])
%     
%     for ii=1:3
%         xlim=get(ax(ii+1),'xlim');
%         ylim=get(ax(ii+1),'ylim');
%         L=plot(ax(ii+1),xlim,[0 0],'-',[0 0],ylim,'-','linew',1.5,'color',[.7 .7 .7]);
%         uistack(L, 'bottom');
%         set(L,'linew',1)
%         set(ax(ii+1),'xlim',xlim,'ylim',ylim)
%     end
%     set(h,'markersize',15)
%     
%     xl=ax(4).XLim;
%     yl=ax(4).YLim;
%     afftyp={'SA1','RA','PC'};
%     for ii=1:3
%         text(xl(1)+diff(xl),yl(1)+(1.2-.12*(ii-1))*diff(yl),...
%             afftyp{ii},'parent',ax(4),'color',affcol(ii),'fontw','bold')
%     end
% end
%% most other parameters
parIdx=[1 8 9 10 13]; l=length(parIdx);
ylabs=parameters.description(parIdx);

am=[1 5 6 7 9];
h=[];
for ii=1:l
    plot(ax(am(ii)),[.6 1.4],nanmean(sap(:,parIdx(ii)))*ones(2,1),'k-',...
        [1.6 2.4],nanmean(rap(:,parIdx(ii)))*ones(2,1),'k-',...
        [2.6 3.4],nanmean(pcp(:,parIdx(ii)))*ones(2,1),'k')
    h(end+1)=plot(ax(am(ii)),0.75+(1:nsa)/2/nsa,sap(:,parIdx(ii)),'.','color',affcol(1));%showIdx(h(end))
    h(end+1)=plot(ax(am(ii)),1.75+(1:nra)/2/nra,rap(:,parIdx(ii)),'.','color',affcol(2));%showIdx(h(end))
    h(end+1)=plot(ax(am(ii)),2.75+(1:npc)/2/npc,pcp(:,parIdx(ii)),'.','color',affcol(3));%showIdx(h(end))
    title(ax(am(ii)),['(' num2str(parIdx(ii)) ')'])
    ylabel(ax(am(ii)),ylabs{ii})
end

set(h,'markersize',12)
text(1,800,'\^\^','par',ax(5),'color',affcol(1),'horiz','center','fontw','bold')
set(ax(am),'xtick',[1 2 3],'xticklabel',{'SA1','RA','PC'},'xlim',[.5 3.5],...
    'xticklabelrot',90)
set(ax(am([1 2 4])),'yscale','log')

set(ax(5),'ylim',[.8 1000],'ytick',[1 10 100 1000])
set(ax(7),'ylim',[1 1000],'ytick',[1 10 100 1000])

%% post spike current coefs
h=[];
h(end+1)=plot(ax(8),-sap(:,11),-sap(:,12),'.','color',affcol(1));
h(end+1)=plot(ax(8),-rap(:,11),-rap(:,12),'.','color',affcol(2));
h(end+1)=plot(ax(8),-pcp(:,11),-pcp(:,12),'.','color',affcol(3));
set(h,'markersize',12)
xlabel(ax(8),'Fast coef. [a.u.]')
ylabel(ax(8),'Slow coef. [a.u.]')
title(ax(8),'(11-12)')
set(ax(8),'xscale','log','yscale','log','box','off',...
    'xlim',[7e-6 20],'ylim',[7e-6 20],...
    'xtick',[.0001 .01 1],'ytick',[.0001 .01 1],'xgrid','on','ygrid','on')

%%
delete(ax(2:4))
%%
clear x y
ll=[2e-2 2; .5e-4 5e2;1e-2  5e3];
t={[.09999 1e0],[.01 1e0 1e2],[1e0 1e3]};
tl={[.1 1e0],[.01 1e0 1e2],[1e0 1e3]};

hs=.25;

tw=.225;
th=.45;

wo=.22;
ho=.50;

b=ho+[.15 .565]*th;
l=wo+[.15 .565]*tw;
h=.385*th; w=.385*tw;
sh=.015*th; sw=.015*tw;
vm=ho+.55*th-sh/2;
hm=wo+.55*tw-sw/2;

for kk=1:3
    
    o=hs*(kk-1);
    ax=[axes('pos',[l(1)+o b(2) w h]),...
        axes('pos',[l(2)+o b(2) w h]),...
        axes('pos',[l(1)+o b(1) w h]),...
        axes('pos',[l(2)+o b(1) w h])];
    ax0=[axes('pos',[hm+o b(2) sw h]),...
        axes('pos',[l(1)+o vm w sh]),...
        axes('pos',[hm+o b(1) sw h]),...
        axes('pos',[l(2)+o vm w sh])];
    set([ax ax0],'nextplot','add')
    signs={@(x,y) x<0&y>0,@(x,y) x>0&y>0,@(x,y) x<0&y<0,@(x,y) x>0&y<0};
    signs0={@(x,y) x==0&y>0,@(x,y) x<0&y==0,@(x,y) x==0&y<0,@(x,y) x>0&y==0};
    
    x{1}=sap(:,kk*2);x{2}=rap(:,kk*2);x{3}=pcp(:,kk*2);
    y{1}=sap(:,kk*2+1);y{2}=rap(:,kk*2+1);y{3}=pcp(:,kk*2+1);
    
    for ii=1:4
        for jj=1:3
            locx=x{jj}; locy=y{jj};
            id=signs{ii}(locx,locy);
            locx(~id)=nan; locy(~id)=nan;
            plot(ax(ii),locx,locy,'.','color',affcol(jj),'markersize',12)
            locx=x{jj}; locy=y{jj};
            id=signs0{ii}(locx,locy);
            locx(~id)=nan; locy(~id)=nan;
            plot(ax0(ii),locx,locy,'.','color',affcol(jj),'markersize',12)
        end
    end
    
    set(ax([1 3]),'yaxisl','right')
    set(ax([3 4]),'xaxisl','top')
    set(ax,'xscale','log','yscale','log','box','off','xgrid','on','ygrid','on')
    set(ax0,'xcolor','none','ycolor','none','xgrid','on','ygrid','on','xtick',0,'ytick',0)
    set(ax([2 4]),'xlim',ll(kk,1:2),'xtick',t{kk},'xticklabels',[])
    set(ax([1 3]),'xlim',-ll(kk,[2 1]),'xtick',-flip(t{kk}),'xticklabels',[])
    set(ax([1 2]),'ylim',ll(kk,1:2),'ytick',t{kk},'yticklabels',[])
    set(ax([3 4]),'ylim',-ll(kk,[2 1]),'ytick',-flip(t{kk}),'yticklabels',[])
    set(ax(1:2),'xticklabels',[])
    set(ax([1 3]),'yticklabels',[])
    set(ax0(1),'ylim',ll(kk,1:2),'yscale','log','xlim',[-1 1])
    set(ax0(3),'ylim',-ll(kk,[2 1]),'yscale','log','xlim',[-1 1])
    set(ax0(2),'xlim',-ll(kk,[2 1]),'xscale','log','ylim',[-1 1])
    set(ax0(4),'xlim',ll(kk,1:2),'xscale','log','ylim',[-1 1])
    
    opt1={'horiz','right','vert','mid','fontsize',7};
    opt2={'horiz','center','vert','mid','fontsize',7};
    for gg=1:length(t{kk})
        text(-5,t{kk}(gg),num2str(tl{kk}(gg)),'par',ax0(1),opt1{:})
        text(-5,-t{kk}(gg),num2str(-tl{kk}(gg)),'par',ax0(3),opt1{:})
        
        text(-t{kk}(gg),10,num2str(-tl{kk}(gg)),'par',ax0(2),opt2{:})
        text(t{kk}(gg),10,num2str(tl{kk}(gg)),'par',ax0(4),opt2{:})
    end
    
    text(0,-10^(log10(ll(kk,1))+log10(ll(kk,2)/ll(kk,1))*1.22),[labels{kk} ' (+)'],'par',ax0(3),'horiz','center','fontsize',11)
    text(-10^(log10(ll(kk,1))+log10(ll(kk,2)/ll(kk,1))*1.22),0,[labels{kk} ' (-)'],'par',ax0(2),'rot',90,'horiz','center','fontsize',11)
    title(ax0(1),titles{kk})
end

%% screen
screenshot(fig,'figs/parameters_plot01_curr','png')
if(screenshot)
    close(fig)
end