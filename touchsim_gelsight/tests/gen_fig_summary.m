% This script is intended to produce a latex file output.  This latex file
% summarize a bunch of figures testing the model:
%  (a) figures generated from the data or a paper figure
%     (might also contain next figure 2 superimposed)
%  (b) original figure generated from the model at the time of writing the script
%  (c) current figure generated from the model currently
%
% You should follow the following formalism:
%  - scriptname.m (script generating the current figure)
%  - scriptname.txt (title [line #1] of the script)
%  - scriptname##_data.eps (a)  (## fig number, e.g. 01,02,03,...)
%  - scriptname##_orig.eps (b)
%  - scriptname##_curr.eps (c)
%  - scriptname##.txt (title [line #1] and comments [line #2] file)
%
% Scripts, text files and figures should be contained in a folder named
% "gen_fig" and which is in the same directory as this script.
function gen_fig_summary(mode)

warning off

if nargin<1
    mode = 'full';
end

t_all = tic;

close all;
clc

% get current version number and machine
[~,rev] = system('hg id');
[~,rev_local] = system('hg id -n');
[~,hostname] = system('hostname');

currp=pwd;
p = mfilename('fullpath');
p = fileparts(p);
cd([p '/gen_fig'])

mlist=dir('*.m');
lm=length(mlist);

fid = fopen('gen_fig_tmp.txt','w');
hasfailed = false;
for ii=1:lm
    % 1. Run the script
    scriptname=mlist(ii).name(1:end-2);
    disp(['Current script : ' num2str(ii) '. ' scriptname])
    if strcmp(mode,'full')
        disp('Running script...')
        try
            t_sub = tic;
            run([p '/gen_fig/' scriptname]);
            close all
            toc(t_sub);
        catch
            warning on
            warning('Script execution failed. Skipping...')
            warning off
            hasfailed = true;
            fprintf(fid,'\\paragraph{Error: execution of %s.m failed.}\\newpage \n\r',strrep(scriptname,'_','\_'));
            continue;
        end
    else
        disp('Skipping execution...')
    end
    
    % 2. Generate latex output with title and description
    txtfid=fopen(sprintf('docs/%s.txt',scriptname),'r');
    if txtfid==-1
        warning('No figure title found.');
        title = '';
    else
        title = fgetl(txtfid);
        fclose(txtfid);
    end
    % 3. loop over each figure
    flist=dir(sprintf('docs/%s*.txt',scriptname));    
    lf=length(flist)-1;
    for jj=1:lf
        disp(['Adding figure ' num2str(jj) ' of ' num2str(lf) '...'])
        txtfid=fopen(sprintf('docs/%s%02d.txt',scriptname,jj),'r');
        if txtfid==-1
            warning('No figure description found.')
            subtitle = '';
            description = '';
        else
            subtitle = fgetl(txtfid);
            description = fgetl(txtfid);
            fclose(txtfid);
        end
        fprintf(fid,'\\paragraph{%s (%s.m)}\\underline{%s} %s\n\r',title,strrep(scriptname,'_','\_'),subtitle,description);
        fprintf(fid,'\\vspace{1cm} \\begin{figure} \\centering \n\r ');
        
        total_figs = double(exist(sprintf('figs/%s%02d_data.eps',scriptname,jj),'file')==2 || exist(sprintf('figs/%s%02d_data.png',scriptname,jj),'file')==2) +...
            double(exist(sprintf('figs/%s%02d_orig.eps',scriptname,jj),'file') == 2) + 1;
        fwid = 1/total_figs;
        if exist(sprintf('figs/%s%02d_data.eps',scriptname,jj),'file')==2 || exist(sprintf('figs/%s%02d_data.png',scriptname,jj),'file')==2
            fprintf(fid,'\\includegraphics[width=%1.1f\\linewidth,height=.6\\textheight,keepaspectratio]{figs/%s%02d_data}',fwid,scriptname,jj);
        else
            disp('No data figure found. Skipping...')
        end
        
        if exist(sprintf('figs/%s%02d_orig.eps',scriptname,jj),'file') == 2
            fprintf(fid,'\\includegraphics[width=%1.1f\\linewidth,height=.6\\textheight,keepaspectratio]{figs/%s%02d_orig}',fwid,scriptname,jj);
        else
            disp('No original model figure found. Skipping...')
            
        end
        
        fprintf(fid,'\\includegraphics[width=%1.1f\\linewidth,height=.6\\textheight,keepaspectratio]{figs/%s%02d_curr}\n\r',fwid,scriptname,jj);
        currplot = dir(sprintf('figs/%s%02d_curr.pdf',scriptname,jj));
        if isempty(currplot)
            currplot = dir(sprintf('figs/%s%02d_curr.png',scriptname,jj));
        end
        fprintf(fid,'\\caption*{Model figure timestamp: %s}\n\r \\end{figure}',currplot.date);
        
        fprintf(fid,'\\newpage\n\r');
    end
end
disp('Tex file generated.')
fclose(fid);

text = fileread('gen_fig_tmp.txt');

fin = fopen('gen_fig_template.tex','r');
fout = fopen('gen_fig.tex','w');
while(~feof(fin))
    s = fgetl(fin);
    if(strfind(s,'#content#'))
        s = strrep(s, '#content#', text);
    elseif(strfind(s,'#datestamp#'))
        s = strrep(s, '#datestamp#', sprintf('\\date{\\today, revision %s (%s on %s)}',rev,rev_local,strtrim(hostname)));
    end
    fprintf(fout,'%s\n\r',s);
end
fclose(fin);
fclose(fout);

% Try compiling to PDF and opening file
timestr=datestr(now,'yyyy_mm_dd_HH_MM_SS');
filenamewithtimestamp=['gen_fig_' timestr '.pdf'];
[status,~] = system('pdflatex -interaction=nonstopmode gen_fig.tex');
if status==0
    disp('PDF file generated.')
    delete('gen_fig.aux','gen_fig.log','gen_fig.tex','gen_fig_tmp.txt')
    movefile('gen_fig.pdf',filenamewithtimestamp)
    if ispc
        system(['START /B .\' filenamewithtimestamp]);
    elseif ismac
        system(['open -a Preview ' filenamewithtimestamp]);
    end
else
    disp('Could not compile to PDF.')
end


cd(currp)

toc(t_all)

if hasfailed
    fprintf(2,'WARNING: Some scripts have failed.\n')
end

warning on
