function y=nansem(x,flag,dim)
%for input, see nanstd

y=nanstd(x,flag,dim)./sqrt(sum(~isnan(x),dim));