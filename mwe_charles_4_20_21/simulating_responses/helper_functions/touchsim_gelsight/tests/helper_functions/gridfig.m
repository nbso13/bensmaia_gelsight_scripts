function ax = gridfig(nr,nc,varargin)
%GRIDFIG Multiple plots drawn in a grid.
%This is an alternative to 'subplot' allowing better control on margins.
%Create an empty fig with multiple subplots.
%Return axes handles, sorted as with the subplot function.
%
% Syntax: ax = gridfig(nr,nc,varargin)
%
% Inputs:
%   nr            number of rows
%   nc            number of columns
%   'parent'      figure handle (if empty creates new fig)
%   'visible'     'on' or 'off', make subaxes visible    (default: on)
%   'shownumber'  T or F, show axis number               (default: F)
%   'spacing'     spacing amplification factor           (default: 1)
%   'spacingv'    vertical spacing (overulled by spacing)(default: 1)
%   'spacingh'    horizont spacing (overulled by spacing)(default: 1)
%   'margin'      bottom-left margin percentage          (default: 0)
%   'marginl'     left margin perc (overulled by margin) (default: 0)
%   'marginb'     bottom margin p  (overulled by margin) (default: 0)
%   'titleline'   T or F, add space for title line       (default: F)
%   'title'       string specifying title                (default: '')
%   'handletype'  'vector' or 'matrix'                   (default: 'vector')
%   'ordering'    'hor' or 'ver' (if handletype='vector')(default: 'hor')
%   'title'       title of the figure                    (default: '')
%   'nextplot'    set axes 'nextplot' property           (default: 'add')
%   'merge'       merge adj subplot, ex:{[1,4], [5,6]}   (default:{})
%      (specifying idx as handletype vector and ordering hor, like subplot)
%      (only available in handletype vector and ordering hor mode)
%      (return full original ax, merged subaxes is refered as the upper
%       left handle, others are deleted and value is set to nan in ax)
%
% Outputs:
%   (none)

p=inputParser;
p.addParamValue('parent',1);
p.addParamValue('visible','on');
p.addParamValue('shownumber',0);
p.addParamValue('spacing',1);
p.addParamValue('spacingv',1);
p.addParamValue('spacingh',1);
p.addParamValue('margin',0);
p.addParamValue('marginl',0);
p.addParamValue('marginb',0);
p.addParamValue('titleline',0);
p.addParamValue('handletype','vector');
p.addParamValue('title','');
p.addParamValue('nextplot','add');
p.addParamValue('ordering','hor');
p.addParamValue('merge',{});

p.parse(varargin{:})
parent=p.Results.parent;
visible=p.Results.visible;
shownumber=p.Results.shownumber;
sp=p.Results.spacing;
spv=p.Results.spacingv;
sph=p.Results.spacingh;
mr=p.Results.margin;
mrl=p.Results.marginl;
mrb=p.Results.marginb;
titleline=p.Results.titleline;
title=p.Results.title;
htype=p.Results.handletype;
nextplot=p.Results.nextplot;
order=p.Results.ordering;
merge=p.Results.merge;

if(~isempty(title) && titleline == 0)
    titleline=1;
end

if(mr~=0)
    mrl=mr;
    mrb=mr;
end

if(sp~=1)
    spv=sp;
    sph=sp;
end

if(nargin<2)
    close all
    parent=figure;
    nr=5;
    nc=5;
end

figure(parent)

if(~isempty(title))
    ha = axes('parent',parent,...
        'Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 1],...
        'Box','off','Visible','off','Units','normalized',...
        'clipping','off','hittest','off');
    text(0.5,.99,title,...
        'HorizontalAlignment','center','VerticalAlignment', 'top',...
        'fontsize',16)
end

ex=.92;
ax=zeros(nr,nc);
for ii=1:nr
    for jj=1:nc
        ax(ii,jj)=axes('par',parent,'unit','norm',...
            'pos',[...
            mrl+(1-mrl)*(jj-(1-.022*sph*nc^ex))/nc ...
            mrb+(1-mrb-.05*titleline)*(ii-(1-.022*spv*nr^ex))/nr...
            (1-mrl)*(1-.028*sph*nc^ex)/nc ...
            (1-mrb-.05*titleline)*(1-.028*spv*nr^ex)/nr ...
            ],...
            'visible',visible);
    end
end

set(ax,'nextplot',nextplot);

if(strcmp(htype,'vector'))
    ax=flipud(ax);
    if(strcmp(order,'hor'))
        ax=ax';
    end
    ax=ax(:);
    
else
    ax=flipud(ax);
end

if(strcmp(htype,'vector') && strcmp(order,'hor'))
    for ii=1:length(merge)
        % do not do anything if length==1
        if(length(merge{ii})==1)
            warning('No merge: single element')
            continue;    
        end
        % put idx into the matrix grid
        mat=false(nc,nr);     mat(merge{ii})=true;
        % test if (define block) and merge
        [x,y]=find(mat);
        if(all(column(mat(min(x):max(x),min(y):max(y)))))
            mat=find(mat);
            axh=ax(mat);
            pos=cell2mat(get(axh,'pos'));
            x=min(pos(:,1));            
            y=min(pos(:,2));
            w=max(pos(:,1))-x+pos(1,3); 
            h=max(pos(:,2))-y+pos(1,4);
            set(axh(1),'pos',[x,y,w,h]);
            delete(axh(2:end));
            ax(mat(2:end))=nan;
        else
            warning('No merge: not forming a block');
        end
    end
end


if(shownumber)
    for ii=1:length(ax)
        text(.5,.5,num2str(ii),'parent',ax(ii))
    end
end

end
