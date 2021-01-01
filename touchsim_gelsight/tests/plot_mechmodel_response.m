clear all; close all; clc

d=logspace(-1,2,6)'; ld=length(d); locs=[d zeros(size(d))];
ap=AfferentPopulation;ap.add_afferents('RA',locs);
f=logspace(0,3,7); lf=length(f);

sinamp=.1;

a=20; % number of cycles
b=100; % number of points per cycle
l=a*b;

dur=1./f*a;
sfreq=f*b;

qnew=zeros(l*lf,ld);
qold=zeros(l*lf,ld);
y=zeros(1,l*lf);

for ii=1:lf
    t=1/sfreq(ii):1/sfreq(ii):dur(ii);
    y(l*(ii-1)+1:l*ii)=sin(2*pi*f(ii)*t)*sinamp+1;

    stimn(ii)=Stimulus(y(l*(ii-1)+1:l*ii)',[0 0],sfreq(ii));
    stimo(ii)=Stimulus(y(l*(ii-1)+1:l*ii)',[0 0],sfreq(ii));
    
    for jj=1:ld
        qnew(l*(ii-1)+1:l*ii,jj)=stimn(ii).propagate(ap.afferents(jj)).trace;
        qold(l*(ii-1)+1:l*ii,jj)=stimo(ii).propagate(ap.afferents(jj)).trace;
    end
end

ampn=zeros(lf,ld);
ampo=zeros(lf,ld);

for ii=1:lf
    idx=l*(ii-1)+1:l*ii;
    idx=idx(7*end/10:9*end/10);
    ampn(ii,:)=(max(qnew(idx,:))-min(qnew(idx,:)))/2/sinamp/1000;
    ampo(ii,:)=(max(qold(idx,:))-min(qold(idx,:)))/2/sinamp;
end


%%
figure('pos',[200 50 1400 940])
set(1,'defaultlinelinewidth',2)
set(1,'DefaultAxesFontSize',14)

subplot(2,2,1)
loglog(f,ampn)
xlabel('Frequency (Hz)')
ylabel('Gain')
title('New model','fontsize',16)
% h=legend(cellstr(num2str(d,'%2.1f')),'location','southeast');
% h = get(h,'Title');
%set(h,'String','Distance','fontsize',14)

subplot(2,2,2)
loglog(d,ampn)
xlabel('Distance (mm)')
ylabel('Gain')
title('New model','fontsize',16)
% h=legend(cellstr(num2str(f','%2.0f')),'location','southwest');
% h = get(h,'Title');
%set(h,'String','Frequency','fontsize',14)

% subplot(2,2,3)
% loglog(f,ampo)
% xlabel('Frequency (Hz)')
% ylabel('Gain')
% title('Old model','fontsize',16)
% h=legend(cellstr(num2str(d,'%2.1f')),'location','southeast');
% h = get(h,'Title');
% set(h,'String','Distance','fontsize',14)
% 
% subplot(2,2,4)
% semilogy(d,ampo)
% xlabel('Distance (mm)')
% ylabel('Gain')
% title('Old model','fontsize',16)
% h=legend(cellstr(num2str(f','%2.0f')),'location','southwest');
% h = get(h,'Title');
% set(h,'String','Frequency','fontsize',14)
