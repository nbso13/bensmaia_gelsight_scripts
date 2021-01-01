ccc

load('D:\Dropbox (INMACOSY)\shares\PeripheralModelTrainingData\training_data.mat')
sf=5000; % sampling freq
%models={@(x) x(:,1:4)    @(x) x(:,1:4)    @(x) x(:,3:6)}; % stimulus
%models={@(x) x(:,1:6)    @(x) x(:,1:6)    @(x) x(:,1:6)}; % stimulus
models={@(x) x(:,1:4)    @(x) x(:,1:4)    @(x) x(:,3:6)}; % stimulus

gw_length=[80 60 40]*sf/1000; % windows length [ms] for computing rate
nfilt=40*sf/1000; % filter length [ms]
tsamp=-(nfilt-1):1:0;

inputs={inputSA,inputRA,inputPC};
cidx=strcmp(aff_class,'SA')+strcmp(aff_class,'RA')*2+strcmp(aff_class,'PC')*3;

% select only noise
% lim=[1.2e5 1.2e5 3.6e5];
% for ii=1:3
%     inputs{ii}=inputs{ii}(1:lim(ii),:);
% end
% for ii=1:17
%     spvec{ii}=spvec{ii}(1:lim(cidx(ii)),:);
%     sp{ii}=sp{ii}(sp{ii}<lim(cidx(ii)),:);
% end

stim=cell(1,3);
scov=cell(1,3);
for ii=1:3
    % SA1, RA and PC stimulus
    stim{ii}=models{ii}(inputs{ii});
    % stim cov
%     X=bsxfun(@(x,y) circshift(x,y,1),permute(stim{ii},[1 3 2]),-tsamp);
%     scov{ii}=mtimesx(permute(X,[2 1 3]),X);
%     for jj=1:size(scov{ii},3)
%         %lambda{ii}(jj)=mean(diag(scov{ii}(:,:,jj)))*5000;
%         [U,S,V] = svd(scov{ii}(:,:,jj));
%         S(S<S(1,1)/20)=0; % low freq
%         S(S>0)=1./S(S>0); % inverse
%         pseudscov{ii}(:,:,jj)=U*S*V';
%     end
end
%% fitting procedure
opts = optimset('Display','off');
paramstrspec=repmat('%.10e ',1,length(tsamp)*size(stim{1},2)+3);

p0=zeros(2,17);
pnlin=zeros(3,17);
classidx=[0 0 0];

disp('parameters = struct;')

% ntest=20;
% rs=zeros(17,ntest,ntest); % rsquare for diff neurons and lambdas
% rss=zeros(17,ntest,ntest); % res sum of square

rs=zeros(17,1); % rsquare for diff neurons and lambdas
rss=zeros(17,1); % res sum of square

for nidx=1:17;
    class=cidx(nidx);
    ndstim=size(stim{class},2);
    t=1:length(stim{class});
    samp=bsxfun(@plus,sp{nidx},tsamp)';
    
    % sta
    stn=interp1(t,stim{class},samp);
    sta=mean(stn,2);
    sta=squeeze(sta);
    
    % 1. STA
    for kk=1:size(sta,2)
        %plin{nidx}(:,kk)=pseudscov{class}(:,:,kk)*sta(:,kk); % this is the STA whitened and regularized
        plin{nidx}(:,kk)=sta(:,kk); % this is the STA 
    end
    
    % 2. lin part
    lprate{nidx}=sum(conv2(stim{class},plin{nidx},'same'),2);
    mrate=conv(spvec{nidx},gw(gw_length(cidx(nidx))),'same');
    [~,edges]=histcounts(lprate{nidx}(:,1));
    
    groups = discretize(lprate{nidx}(:,1),edges);
    grp_mrate=accumarray(groups, mrate, [], @mean);
    grp_lprate=accumarray(groups, lprate{nidx}, [], @mean);
   
    % 3. fit non-lin part
    p0(:,nidx)=regress(mrate,[ones(size(lprate{nidx},1),1) lprate{nidx}]);
    pnlin(:,nidx)=lsqcurvefit(@GLM_nonlinfun,[0 p0(2,nidx) 1],lprate{nidx},mrate,[-inf -inf 0],[inf inf 5],opts);
    %pnlin(:,nidx)=lsqcurvefit(@GLM_nonlinfun,[0 1 1],grp_lprate,grp_mrate,[-inf -inf 0],[inf inf 5],opts);
    
    nlprate{nidx}=GLM_nonlinfun(pnlin(:,nidx),lprate{nidx});
    grp_nlprate=accumarray(groups, nlprate{nidx}, [], @mean);
    
    rs(nidx)=rsquare(mrate,nlprate{nidx});
    rss(nidx)=nansum((mrate-nlprate{nidx}).^2);
        if 1
        figure(2)
        subplot(3,6,nidx)
        maxval=max([grp_nlprate;grp_mrate]);
        plot([0 maxval],[0 maxval],'k--',...
            grp_nlprate,grp_mrate,'.')
        title(aff_class{nidx})
    end
    classidx(cidx(nidx))=classidx(cidx(nidx))+1;
    fprintf(['parameters.' lower(aff_class{nidx}) '{' num2str(classidx(cidx(nidx))) '} = ['...
        paramstrspec '];\n'],[plin{nidx}(:);pnlin(:,nidx)])
    %         end
    %     end
    
    spvec_pred{nidx}=nlprate{nidx}>rand(size(nlprate{nidx}));
    
    
    % test if I recover STA
    if 0
        figure(3)
        subplot(3,6,nidx)
        plot(tsamp,plin{nidx}); hold on
        set(gca,'colororderindex',1)
        sta2=zeros([size(sta) 20]);
        for tt=1:40
            
            spikes=find(spvec_pred{nidx});
            samp=bsxfun(@plus,spikes,tsamp)';
            sta2(:,:,tt)=squeeze(nanmean(interp1(t,stim{class},samp),2));
        end
        for kk=1:size(sta,2)
            plotci(tsamp,pseudscov{class}(:,:,kk)*squeeze(sta2(:,kk,:))); 
        end
        hold off
    end
end
return
%%
for ii=1:17
    subplot(3,6,ii)
    plot(cumsum(spvec{ii}),'k'); hold on
    pred=zeros(length(spvec{ii}),10);
    for jj=1:10
        pred(:,jj)=(fr_nlinpred{ii})>rand(size(fr_nlinpred{ii}));
    end
    plot(cumsum(pred),'col',[.7 .7 .7])
end
%%
rs(rs<0)=0;
for ii=1:17
    subplot(4,5,ii)
    loc=squeeze(rs(ii,:,:));
    imagesc(loc); hold on
    [val1,idx1]=max(loc,[],1);
    [val2,idx2]=max(val1,[],2);
    plot(idx2,idx1(idx2),'rx')
    title([aff_class{ii} ' ' num2str(val2)])
    lambda(ii,:)=2.^([idx1(idx2),idx1(idx2),idx2,idx2]*2);
    colorbar
end

class={'SA','RA','PC'};
for ii=1:3
    subplot(4,5,17+ii)
    loc=squeeze(mean(rs(strcmp(aff_class,class{ii}),1:10,1:10),1));
    imagesc(loc); hold on
    [val1,idx1]=max(loc,[],1);
    [val2,idx2]=max(val1,[],2);
    plot(idx2,idx1(idx2),'rx')
end

%% TEST
