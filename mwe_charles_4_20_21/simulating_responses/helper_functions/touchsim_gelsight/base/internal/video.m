function video(r,slowdown,pop,filename)
% video(r,slowdown,pop,filename)

if nargin<4
    filename = ['v' strrep(num2str(fix(clock())),' ','')];
end

if nargin<3
    pop = 'sep';
end

if nargin<2
    slowdown = 1;
end

% 30 frames/s -> 33.3 ms time windows
bin_time=33.3/slowdown;

r_time = r.psth(bin_time);
r_time(r_time==0) = NaN;

tt = linspace(0,r.stimulus.duration,size(r_time,1));

r_timeSA = r_time;
r_timeSA(:,r.affpop.iRA | r.affpop.iPC) = NaN;
r_timeRA = r_time;
r_timeRA(:,r.affpop.iSA1 | r.affpop.iPC) = NaN;
r_timePC = r_time;
r_timePC(:,r.affpop.iRA | r.affpop.iSA1) = NaN;

figure(1)
set(1,'color','w');
if ~strcmp(pop,'all')
    set(1,'pos',[0 0 1200 500])
end

v = VideoWriter([filename '.avi']);
open(v);

audio_sf=30*1470;
t=1/audio_sf:1/audio_sf:bin_time/1000;
spkshape=[1-exp(-t(1:50)/.0005) exp(-t(51:end)/.002)];
sasound=sin(2*pi*t*2000).*spkshape;
rasound=sin(2*pi*t*4000).*spkshape;
pcsound=sin(2*pi*t*8000).*spkshape;

anySA=nansum(r_timeSA,2)>0;
anyRA=nansum(r_timeRA,2)>0;
anyPC=nansum(r_timePC,2)>0;


v2=vision.VideoFileWriter('Filename',[filename '_with_audio.avi'],...
    'FrameRate',30,'AudioInputPort',1);
% v2.AudioCompressor='MJPEG Compressor';
v2.VideoCompressor='MJPEG Compressor';

for t=1:size(r_time,1)
    clf
    
    if strcmp(pop,'all')
        col = repmat([.8 .8 .8],r.affpop.num,1) - repmat(r_time(t,:)',1,3).*repmat(([.8 .8 .8]-affcol(3))./max(r_time(:)),r.affpop.num,1);
        plot_hand(gca,'names',false,'axes',false,'centers',false,'afferents',r.affpop,'color',col);
        title(num2str(t))
    else
        subplot(3,3,[1,4])
        col = repmat([.8 .8 .8],r.affpop.num,1) - repmat(r_timeSA(t,:)',1,3).*repmat(([.8 .8 .8]-affcol(1))./max(r_timeSA(:)),r.affpop.num,1);
        plot_hand(gca,'names',false,'axes',false,'centers',false,'afferents',r.affpop,'color',col);
        title('SA1')
        
        subplot(3,3,[2,5])
        col = repmat([.8 .8 .8],r.affpop.num,1) - repmat(r_timeRA(t,:)',1,3).*repmat(([.8 .8 .8]-affcol(2))./max(r_timeRA(:)),r.affpop.num,1);
        plot_hand(gca,'names',false,'axes',false,'centers',false,'afferents',r.affpop,'color',col);
        title('RA')
        
        subplot(3,3,[3,6])
        col = repmat([.8 .8 .8],r.affpop.num,1) - repmat(r_timePC(t,:)',1,3).*repmat(([.8 .8 .8]-affcol(3))./max(r_timePC(:)),r.affpop.num,1);
        plot_hand(gca,'names',false,'axes',false,'centers',false,'afferents',r.affpop,'color',col);
        title('PC')
        
        subplot(3,3,[7,8,9])
        plot(linspace(0,r.stimulus.duration,size(r.stimulus.trace,1)),r.stimulus.trace,'LineWidth',1.5)
        line([tt(t) tt(t)],get(gca,'ylim'),'Color','k','LineWidth',2)
        xlim([0 r.stimulus.duration])
        box off
        %axis off
        xlabel('Time [s]')
        ylabel('Indentation [mm]')
    end
    
    % add frame to AVI
    F = getframe(1);
    writeVideo(v,F)
    
    audio=anySA(t)*sasound+anyRA(t)*rasound+anyPC(t)*pcsound;
    step(v2,F.cdata,[audio' audio'])

    % add frame to GIF
    im = frame2im(F);
    [imind,cm] = rgb2ind(im,256);
    if t == 1;
        imwrite(imind,cm,[filename '.gif'],'gif', 'Loopcount',inf,'DelayTime',0.0333);
    else
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','DelayTime',0.0333);
    end
end
close(v)
release(v2)
close(1)
