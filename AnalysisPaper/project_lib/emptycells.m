% c = emptycells(c)
% 
% Konvertiert leere Zellen zu leeren Strings
%  z.B. c{1}=[] wird zu c{1}=''
function c = emptycells(c);
for i=1:length(c)
    if(isempty(c{i}))
        c{i}='';
    end
end
