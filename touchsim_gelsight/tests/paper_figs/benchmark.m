ccc
pins=[1 3 10 30 100 300 1000];
recs=[1 3 10 30 100 300 1000];
sf=[100 300 1000 3000 10000];
[x,y]=meshgrid(-5:.1:5);
lp=length(pins);
lr=length(recs);
lsf=length(sf);

t=zeros(lp,lr,4);
tottime_summary=zeros(lp,lr,lsf);
lab={'stim','affpop','resp','all'};

for kk=1:lsf
    for ii=1:lp
        for jj=1:lr
            npin=pins(ii);
            nrec=recs(jj);
            trace=randn(sf(kk),npin);
            
            timall=tic;
            s=Stimulus(trace,[x(1:npin)' y(1:npin)'],sf(kk));
            t(ii,jj,1)=toc(timall);
            
            timaffpop=tic;
            ap=AfferentPopulation;
            ap.add_afferents('RA',[x(1:nrec)' y(1:nrec)']);
            t(ii,jj,2)=toc(timaffpop);
            
            
            timresp=tic;
            rc=ap.response(s);
            t(ii,jj,3)=toc(timresp);
            t(ii,jj,4)=toc(timall);
        end
    end
%     summary=[0 recs;pins' t(:,:,4)];
%     summary=reshape(cellstr(num2str(summary(:),3)),size(summary));
%     summary{1}='\#pins\textbackslash\#recs';
%     cell2latex(['../doc/model_perf/perf_' num2str(sf(kk)) 'hz.tex'],summary,1,'lightgray',num2str(sf(kk)))
    tottime_summary(:,:,kk)=t(:,:,4);
end

%%
close all
figure(1); set(1,'pos',[680   224   560   754])
ax=zeros(1,5);
for ii=1:lsf
    ax(ii)=subplot(3,2,ii);
    imagesc(1:lp,1:lr,tottime_summary(:,:,ii),'parent',ax(ii))
    title(sprintf('sf: %d Hz',sf(ii)))
    xlabel('# pins')
    ylabel('# receptors')
    box off
end
set(ax,'clim',[0 1],'xtick',1:lp,'xticklabels',pins,...
    'ytick',1:lr,'yticklabels',recs)
subplot(3,2,6)
set(gca,'visible','off')
c=colorbar; c.Label.String = 'Duration [s]';

screenshot('benchmark')