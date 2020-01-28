% h = abplot(a, b, varargin)
%   Plottet eine Line y=ax+b
%   nimmt get(gca,'YLim') und get(gca,'YLim') als x- und y-Werte
% 
%   a       Steigung
%   b       y-Achsenabschnitt
%   varargin wird plot übergeben
%           Default: 'k' 
% 
%   h       handle des plot-Befehls

function h = abplot(a,b,varargin)
if(nargin==0)
    disp('h = abplot(a,b,varargin)')
    return
end
if(nargin==1 | isempty(b))
    b = a(2);
    a = a(1);
end

xl = get(gca,'XLim');
% ylim = get(gca,'YLim');

if(nargin<3)
    args{1} = 'k';
else
    args = varargin;
end

if(strcmp(get(gca,'NextPlot'),'replace'))
    holdon = 1;
else
    holdon = 0;
end

if(strcmp(get(gca,'XScale'),'log')==1)
    x = logspace(log10(xl(1)),log10(xl(2)),100);
else
    x = linspace(xl(1),xl(2),100);
end
if(isinf(a))
    y = get(gca,'YLim');
    x = [b,b];
else
    y = a*x+b;
end

if(holdon==1)
    hold on
    h=plot(x,y,args{:});
    hold off
else
    h=plot(x,y,args{:});
end    


