% clean up html directory
delete docs/html/*.html
delete docs/html/*.png

% publish all files in docs directory
mlist=dir('docs/*.m');
lm=length(mlist);
for i=1:lm
    publish(mlist(i).name);
    close all
end
