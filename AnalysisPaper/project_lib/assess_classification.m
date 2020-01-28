% ass = assess_classification(x,y)
% 
function ass = assess_classification(x,y,fast)
if nargin<3
    fast = false;
end

if ~fast

    ass.PPV = sum(x==1 & y==1) / sum(y==1);
    
    ass.TPR = ass.PPV;%sum(x==1 & y==1)/(sum(y==1));
    ass.FPR = sum(x~=1 & y==1)/(sum(x~=1));
    
    ass.TNR = sum(x~=1 & y~=1)/(sum(y~=1));
    ass.FNR = sum(x==1 & y~=1)/(sum(y==1));
    
    ass.FDR = sum(x~=1 & y==1) / sum(y==1);
    ass.ACC = sum(x==y) / length(x);

    ass.Sens = sum(x==1 & y==1) / sum(x==1);
    ass.Spec = sum(y~=1 & x~=1) / sum(x~=1);
    ass.recall = ass.Sens;
    ass.F1 = 2*sum(x==1 & y==1)/(2*sum(x==1 & y==1)+sum(x~=1 & y==1) + sum(y~=1 & x==1));

else
    ass.F1 = 2*sum(x==1 & y==1)/(2*sum(x==1 & y==1)+sum(x~=1 & y==1) + sum(y~=1 & x==1));
end

