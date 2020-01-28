% Ermittelt die Levels des discreten Faktors
% [l,val] = levels(df);
% 
%   l   levels
%   val Werte: Welchem Level entspricht der Vektor?
% 

function [l,val] = levels(df);

% val = NaN*ones(1,length(df));
if(iscell(df) & ~isempty(df))
	l{1} = df{1};
	for i=2:length(df)
        if(isempty(find(strcmp(df{i},l))))
            l{length(l)+1} = df{i};
        end
	end
    
    for i=1:length(l)
        val(find(strcmp(l{i},df))) = i;
    end
elseif(~isempty(df))
	l(1) = df(1);
	for i=2:length(df)
        if(isempty(find(df(i)==l)))
            l(length(l)+1) = df(i);
        end
	end
    l = sort(l);
        
    nan = find(isnan(l));
    if(length(nan)>1)
%         val = val(1:(end-length(nan)+1));
        l = l(1:(end-length(nan)+1));
    end
    for i=1:length(l)
       if(~isnan(l(i))) 
           val(i) = nansum(df==l(i));
       else
           val(i) = sum(isnan(df));
       end
    end
elseif(iscell(df)) % leer und Zelle
    l = {};
    val = [];
else % leer
    l = [];
    val = [];
end    

