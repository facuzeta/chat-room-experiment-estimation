function [r,p]=nancorr(x,y)

iout = isnan(x) | isnan(y) | abs(x)==Inf | abs(y)==Inf;
x(iout)=[];y(iout)=[];
% [r,p]=corr(x,y,'type','spearman');
[r,p]=corr(x,y);

end