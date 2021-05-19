a = AfferentStream('RA');
tic
for i=1:1000
    r = a.response(1,30);
end
t = toc;
fprintf(['Update time: ' num2str(t) ' ms.\n'])
