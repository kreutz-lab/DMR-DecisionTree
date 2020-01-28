% Filtert die Gene und gibt nur noch diejenigen Gene zurück, die das
% Filterkriterium 'filter' erfuellen.
%   filter  boolscher Vektor oder Indizes

function Efilter = FilterGene(E,filter);
ngene = get(E,'ngene');
if( sum(filter==0  | filter==1) == length(filter))
    if(length(filter) == get(E,'ngene'))
        filter = find(filter);  % wandelt boolsche Varialble in Indizes um.
    elseif(length(filter)>1)
        error('FilterGene.m: length(filter) must be equal number of genes.')
    end
end

Efilter = E;


ma      = get(E,'arrays');
for i=1:length(ma)
    data  = get(ma{i},'data');
    ma{i} = set(ma{i},'data',data(filter));
end
Efilter = set(Efilter,'arrays',ma);

IDs     = get(E,'IDs');
Efilter = set(Efilter,'IDs',IDs(filter));

unigeneID = get(E,'unigeneID');
if(~isempty(unigeneID))
    Efilter   = set(Efilter,'unigeneID',unigeneID(filter));
end
namen     = get(E,'namen');
if(~isempty(namen))
    Efilter   = set(Efilter,'namen',namen(filter));
end
symbol    = get(E,'symbol');
if(~isempty(symbol))
    Efilter   = set(Efilter,'symbol',symbol(filter));
end
llid      = get(E,'llid');
if(~isempty(llid))
    Efilter   = set(Efilter,'llid',llid(filter));
end
chromosome= get(E,'chromosome');
if(~isempty(chromosome))
    Efilter   = set(Efilter,'chromosome',chromosome(filter));
end
chrloc= get(E,'chromloc');

if(~isempty(chrloc))
    Efilter   = set(Efilter,'chromloc',chrloc(filter,:));
end
go        = get(E,'go');
if(~isempty(go))
    Efilter   = set(Efilter,'go',go(filter));
end
goids        = get(E,'goids');
if(~isempty(goids))
    Efilter   = set(Efilter,'goids',goids(filter));
end
gonames        = get(E,'gonames');
if(~isempty(gonames))
    Efilter   = set(Efilter,'gonames',gonames(filter));
end
pres        = get(E,'present');
if(~isempty(pres))
    Efilter   = set(Efilter,'present',pres(filter,:));
end

goids        = get(E,'goids');
if(~isempty(goids))
    Efilter   = set(Efilter,'goids',goids(filter,:));
end

go        = get(E,'go');
if(~isempty(go))
    Efilter   = set(Efilter,'go',go(filter,:));
end

try
    cn        = get(E,'cn');
    if(~isempty(cn))
        Efilter   = set(Efilter,'cn',cn(filter));
    end
catch
    warning(lasterr)
end

if(isfield(E.annotation,''))
    if(~isempty(E.annotation.enzyme))
        Efilter = set(Efilter,E.annotation.enzyme(filter));   
    end
end

enzyme = get(E, 'enzyme');
if(~isempty(enzyme))
    Efilter = set(Efilter,'enzyme',enzyme(filter));
end

pathways = get(E, 'pathways');
if(~isempty(pathways))
    Efilter = set(Efilter,'pathways',pathways(filter));
end

try
    cont = get(E,'container');
    fn = fieldnames(cont);
    for i=1:length(fn)
        if(size(cont.(fn{i}),1)==ngene)
            cont.(fn{i}) = cont.(fn{i})(filter,:);
        end
    end
    Efilter = set(Efilter,'container',cont);
catch
    warning(lasterr)
end