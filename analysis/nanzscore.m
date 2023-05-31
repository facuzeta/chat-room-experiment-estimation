function z=nanzscore(y,type)
if nargin<2
    type='norm';
end
   switch type
       case 'norm'
           z=(y-nanmean(y))./nanstd(y);
       case 'nonpar'
           z=(y-nanmedian(y))./mad(y,1);
   end
end