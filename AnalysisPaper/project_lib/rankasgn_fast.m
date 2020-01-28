% [r,rdown,rup] = rankasgn_fast(x)
% 
%   Ist schneller als rankasgn wenn, selten identische Werte vorkommen.
%   Unentschieden führt zu gemittelten Raengen, aber:
% 
%   rdown        Bei Unentschieden der minimale Rang davon
%   rup         Bei Unentschieden der maximale Rang davon

function [r,rdown,rup] = rankasgn_fast(x)

s = size(x);
lin = prod(s);
if(min(s)>1)
    error('rankasgn_fast.m: arrays required.')
end

[insort,rf] = sort(x);
if(nanmin(diff(insort))>0)   
    r(rf) = 1:lin;
    rdown = r;
    rup = r;
else
    r = Rankasgn(x);
    if(nargout>1)
        rdown = Rankasgn(x,@min);
        rup = Rankasgn(x,@max);
    end
end
if(size(x,1)==1 & size(r,1)>1)
    r = r';
end
if(size(x,2)==1 & size(r,2)>1)
    r = r';
end

% inf:
ind = find(x==inf);
if(numel(ind)>1)
    r(ind) = mean(r(ind));
    if(nargout>1)
        rdown(ind) = min(r(ind));
        rup(ind) = max(r(ind));
    end
end
% -inf:
ind = find(x==-inf);
if(numel(ind)>1)
    r(ind) = mean(r(ind));
    if(nargout>1)
        rdown(ind) = min(r(ind));
        rup(ind) = max(r(ind));
    end
end
% NaN
ind = find(isnan(x));
if(numel(ind)>1)
    r(ind) = mean(r(ind));
    if(nargout>1)
        rdown(ind) = min(r(ind));
        rup(ind) = max(r(ind));
    end
end
