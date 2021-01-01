function surround_suppression

% set up afferent population
aSA = affpop_single_models([0 0],'SA1');
aRA = affpop_single_models([0 0],'RA');

% set up stimulus and generate responses
cpin = [0 0];
pins = [ 1 -0.5;
         1  0.5;
         0 -1;
         0  1;
        -1 -0.5;
        -1  0.5];

s{1} = stim_ramp(0.1,0.2,cpin,[],0.0125,[],0.25,1.5);
rSA(:,1) = aSA.response(s{1}).rate;
rRA(:,1) = aRA.response(s{1}).rate;
for n=1:6
    perms = nchoosek(1:6,n);
    rSA_tmp = [];
    rRA_tmp = [];
    for i=1:size(perms,1)
        s{n+1}(i) = stim_indent_shape([cpin; pins(perms(i,:),:)],stim_ramp(0.1,0.2,cpin,[],0.0125,[],0.25,1.5));
        rSA_tmp(:,i) = aSA.response(s{n+1}(i)).rate;
        rRA_tmp(:,i) = aRA.response(s{n+1}(i)).rate;
    end
    rSA(:,n+1) = mean(rSA_tmp,2);
    rRA(:,n+1) = mean(rRA_tmp,2);
end

%% plots results

close all

figure(1)
set(1,'pos',[100 100 800 400])

subplot(121)
semilogy(1:7,mean(rSA),'-ok','LineWidth',2)
%ylim([1 40])
%set(gca,'ytick',[4 10 40])
box off
xlabel('Number of probes')
title('SA1')

subplot(122)
semilogy(1:7,mean(rRA),'-ok','LineWidth',2)
%ylim([1 40])
%set(gca,'ytick',[4 10 40])
box off
xlabel('Number of probes')
title('RA')

screenshot(1,'figs/surround_suppression01_curr','pdf')
close(1)
