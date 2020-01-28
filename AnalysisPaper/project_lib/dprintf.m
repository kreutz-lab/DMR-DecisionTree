% out = dprintf(in)
% 
%  Konvertiert Zahlen, so dass die Ausgabe gut ist.

function out = dprintf(in)
if(iscell(in))
    out = cell(1,length(in));
    for i=1:length(in)
        out{i} = dprintf(in{i});
    end
elseif(length(in)>1)
    out = cell(size(in));
    for i=1:length(in)
        out{i} = dprintf(in(i));
    end    
elseif(~isempty(in))
    if(ischar(in))
        out = sprintf('%s',in);
    elseif(isinf(in))
        if(in<0)
            out = '-inf';
        else
            out = 'inf';
        end
    elseif(isnumeric(in) & in>10)
        out = sprintf('%1.0f',in);
    elseif(isnumeric(in))
        out = sprintf('%1.2g',in);
    else
        in
        error('dprintf.m: in hat falsches Format.')
    end
    out = strrep(out,'e00','e');
    out = strrep(out,'e0','e');
    out = strrep(out,'e-00','e-');
    out = strrep(out,'e-0','e-');
    out = strrep(out,'e+00','e');
    out = strrep(out,'e+0','e');
else
    out = in;
end
