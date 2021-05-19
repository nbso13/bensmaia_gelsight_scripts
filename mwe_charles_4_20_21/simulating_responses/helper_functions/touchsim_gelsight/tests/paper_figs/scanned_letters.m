function scanned_letters

%% create letters
fprintf('create letter...'); tic

letsiz=5; % in mm
pins_per_mm=10;
speed=50; % mm/s
sf=500;%Hz;
stepsiz=speed/sf;
[x,y] = meshgrid((-200:200)./400.*letsiz);
[X,Y] = meshgrid(-letsiz/2:1/pins_per_mm:letsiz/2);

letters='ABCDEFGHIJKLM';

% figure to generate letters layout
h = figure('position',[680   558   401   401]);
set(h,'Visible','off')
axis off
allLetters=zeros(letsiz*pins_per_mm+1);
for ll=1:length(letters)
    text(0.5,0.55,letters(ll),'FontSize',400,'horiz','center','vertic','middle')
    set(gca,'pos',[0 0 1 1])
    A = double(rgb2gray(print2array(h)));
    cla
    
    A = interp2(x,y,A,X,Y);
    A=A==0;
    A(:,sum(A)==0)=[];
    allLetters=[allLetters A zeros(letsiz*pins_per_mm+1,30)];
end
allLetters=[allLetters zeros(letsiz*pins_per_mm+1)];
%allLetters=imfilter(allLetters,fspecial('average',2),'replicate');
close(h)
toc
%% create trace
fprintf('create trace and stimulus...'); tic
[allX,allY]=meshgrid(-letsiz/2:1/pins_per_mm:920/pins_per_mm,...
    -letsiz/2:1/pins_per_mm:letsiz/2);
F = griddedInterpolant(allX',allY',allLetters','linear');
tmax=size(allLetters,2)-letsiz*pins_per_mm;
trace=zeros(tmax,size(allLetters,1)^2);

nsteps=(max(allX(:))-min(allX(:)))/stepsiz;
xy=[Y(:),X(:)]; R=hypot(X(:),Y(:));
for ii=1:nsteps
    lett=F(X'+ii*stepsiz,Y')';
    trace(ii,:)=real(sqrt(2.5^2-R.^2)).*lett(:)...
        +abs(randn(size(R)))*.4; % some 'motor' noise
end
out=find(sum(trace,1)==0);

xy(out,:)=[];
trace(:,out)=[];

stim=Stimulus(trace,xy,sf,1/pins_per_mm/2);
toc

%%
fprintf('compute response...'); tic
affpos=-letsiz/2:.1:letsiz/2;
ap=AfferentPopulation;
for ii=1:length(affpos)
    aploc=affpop_single_models([affpos(ii),0]);
    ap.afferents=[ap.afferents aploc.afferents];
end
rc=ap.response(stim);
spikes=cell(length(affpos),17);
for ii=1:17
    spikes(:,ii)={rc.responses(ii:17:end).spikes};
end
spikes=flipud(spikes);

% clip
spikes=cellfun(@(x) x(x>.014)-.014,spikes,'uni',0);
spikes=cellfun(@(x) x(x<1.7),spikes,'uni',0);
toc

%%
figure(1)
set(1,'pos',[100 250 860 320])
% profile
ax(4)=subplot(411);
imagesc(allLetters)

% SA1
ax(1)=subplot(412);
plot_spikes(spikes(:,2),'linew',.5); 
xlabel(''); ylabel('SA1')

% RA (aligned for delay)
ax(2)=subplot(413);
plot_spikes(cellfun(@(x) x-.007,spikes(:,5),'uni',0),'linew',.5); 
xlabel(''); ylabel('RA')

% PC (aligned for delay)
ax(3)=subplot(414);
plot_spikes(cellfun(@(x) x-.005,spikes(:,14),'uni',0),'linew',.5); 
xlabel('Position [mm]'); ylabel('PC')

set(ax(1:3),'xminortick','on','xtick',0:.2:1.7,'xticklabel',0:10:90,...
    'xlim',[-.01 1.7],'ylim',[-10 60],'ytick',[0:50:50],...
    'yticklabel',0:5:5,'yminortick','on')
set(ax(1:2),'xticklabel',[])
set(ax(4),'visible','off','xlim',[30 890],'ylim',[-7 57])
colormap(flipud(gray))
for ii=1:4
    set(ax(ii),'pos',get(ax(ii),'pos')+[-.07 .045 .15 0.05])
end

screenshot scanned_letters pdf