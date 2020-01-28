% r = rankasgn(x,tiefun)
%
% Assigns ranks to a single vector of values.
%
%     tiefun    function for ranks of ties
%               Default: @mean
%               Here, their mean ranks (=midranks) are assigned.
%
%     Usage: r = rankasgn(x)
%
%           x = [n x 1] data vector.
%           ----------------------------------
%           r = corresponding vector of ranks.
%

function r = rankasgn(x,tiefun)
if(~exist('tiefun','var') | isempty(tiefun))
    tiefun = @mean;
end

[n,c] = size(x);

rowvect = 0;
if (n==1 & c>1)
    x = x';
    [n,c] = size(x);
    rowvect = 1;
end;

if (c>1)
    error('RANKASGN: input vector required');
end;

r = zeros(n,1);               % Allocate output vector

[x,index] = sort(x);          % Sort and get indices to x
i = 1;
d = diff(x)~=0; % Clemi
while (i<n)
    %     if (x(i+1) ~= x(i))           % Not a tie
    if (d(i))           % Clemi
        r(index(i)) = i; % Stash rank
        i = i+1;
    else                          % Tie
        J = n; % is overwritten in case x(i)~=x(j)
        for j = (i+1):n               % How far does it go?
            if (x(i) ~= x(j))
                J = j-1;
                break
            end
        end
        %       midrank = mean(i:j);
        %       for k = i:j                   % Enter midrank into all the tied entries
        %         r(index(k)) = midrank;
        %       end;
        r(index(i:J)) = feval(tiefun,i:J);  %Clemi
        i = J+1;
    end
end
if (r(index(n)) == 0)           % Check last entry
    r(index(n)) = n;
end;

if (rowvect)
    r = r';
end;

if(sum(isnan(x))>0)
    r(index(isnan(x))) = feval(tiefun,r(index(isnan(x)))); %Clemens.
end


