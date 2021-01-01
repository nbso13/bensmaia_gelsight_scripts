function perf=LDAclassL1O(featmatrix,grps,reps)

%check variances
idvar=var(featmatrix)>1e-10; 
%sum(idvar)/length(idvar)

if(~any(idvar))
    perf=ones(1,max(reps))/length(unique(grps));
else
    featmatrix=featmatrix(:,idvar);
    perf=zeros(1,max(reps));
    for ii=1:max(reps)
        train=reps~=ii;
        
        % PCA pre-LDA (regularization for no-cariance within class issues)
        coeff = pca(featmatrix(train,:));
        trainmat=featmatrix(train,:)*coeff;
        try
            obj = fitcdiscr(trainmat,grps(train),...
                'discrimtype','Lin');
        catch
            obj = fitcdiscr(trainmat,grps(train),...
                'discrimtype','pseudoLin');
            disp('used pseudoLinear')
        end
        [~,score] = predict(obj,featmatrix(~train,:)*coeff);
        perf(ii)=mean(diag(score));
    end
end