% function screenshot(filenamewithoutextension)
% OR
% function screenshot(fig,filenamewithoutextension)
% OR
% function screenshot(fig,filenamewithoutextension,format)
% OR
% function screenshot
% OR
% function screenshot(flag)
%
% Takes a screenshot of the current fig (OR figure 'fig')
% - specify 'filenamewithoutextension' without extension
% - you can specify a figure handles with 'fig'
% - you can choose either 'pdf' or 'png' as an output format (default: png)
% - output dimension are the current figure dimensions
%
% Using flag parameter, you can activate/disactivate screenshots
%  flag=0 => disables screenshots
%  flag=1 => enables screenshots for 5 min (and then sets flag to 0)
%  flag=2 => screenshot are always enabled (default behavior)
%  flag=3 => switch status from 0 to 1 or 1 to 0
% Flag status remains even after closing matlab

function status=screenshot(flagorfig,filnamwithoutext,format)

% setup/retrieve global variables
global glob_screenshotflag
global glob_screenshottim
if(isempty(glob_screenshotflag))
    glob_screenshotflag=2;
end

% One numeric argument: modify screenshot flag
if(nargin==1 && ~ischar(flagorfig))
    % first delete all prev callback timers
    if(isobject(glob_screenshottim))
        stop(glob_screenshottim)
        delete(glob_screenshottim)
        glob_screenshottim=[];
    end
    if(flagorfig==3)
        if(glob_screenshotflag>0)
            flagorfig=0;
        else
            flagorfig=1;
        end
    end
    if(flagorfig==0)
        disp('Screenshots disabled')
    end
    if(flagorfig==1)
        % create a timer to disactivate screenshots later
        tim = timer;
        tim.ExecutionMode='singleShot';
        tim.TimerFcn = @(x,y) screenshot(0);
        tim.StartDelay = 300; % (5 min)
        disp('Screenshots have been activated for 5 minutes')
        glob_screenshottim=tim;
        start(tim)
    end
    if(flagorfig==2)
        disp('Screenshots activated all time.')
    end
    glob_screenshotflag=flagorfig;
    return
end

% no argument: get status
if(nargin==0)
    status=glob_screenshotflag;
    return
end

% Other cases, check if going further
if(glob_screenshotflag==0)
    disp('no screenshot')
    return
end

% default behavior: png
if(nargin<3)
    format='png';
end

% if figure handle not provided, shift arguments
if(~ishandle(flagorfig))
    if(nargin<2)
        format='png';
    else
        format=filnamwithoutext;
    end
    filnamwithoutext=flagorfig;
    flagorfig=gcf;
end

% hide all uicontrols (in case of GUI's)
uic=findobj(flagorfig,'type','UIControl');
set(uic,'visible','off')

% save figure properties
figprop=get(flagorfig);

% set figure properties for screenshot
set(flagorfig,'PaperPositionMode','auto','Units','inches','PaperUnits','inches');
pos=get(flagorfig,'pos');
set(flagorfig,'PaperSize',[pos(3), pos(4)])
set(flagorfig,'color','w') % set white background
set(flagorfig,'InvertHardcopy','off') % keep all other visual porperties

% print command depending on output format
switch format
    case 'pdf', print(flagorfig,[filnamwithoutext '.pdf'],'-dpdf','-painters')
    case 'png', print(flagorfig,[filnamwithoutext '.png'],'-dpng','-r300')
    otherwise, disp('unknown format, abord...');
end

% show uicontrols back
set(uic,'visible','on');

% reset figure properties to normal (after removing read-only fields)
set(flagorfig,'PaperPositionMode',figprop.PaperPositionMode)
set(flagorfig,'Units',figprop.Units)
set(flagorfig,'PaperSize',figprop.PaperSize)
set(flagorfig,'color',figprop.Color)
set(flagorfig,'InvertHardcopy',figprop.InvertHardcopy)

disp(['SCREENSHOT SAVED : ' filnamwithoutext])

end