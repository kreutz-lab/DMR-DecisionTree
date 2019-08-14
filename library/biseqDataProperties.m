function data = DataProperties(data)

data.isF = zeros(size(data.strand));
data.isF(strmatch('F',data.strand,'exact')) = 1;
data.isF(strmatch('+',data.strand,'exact')) = 1;
data.isF(strmatch('C',data.strand,'exact')) = 1;  % Aug 2018: das nehme ich für Arabidopsis einfach mal an, da zeigt "C" und "G" den strand an, für die Analyse spielts keine Rolle
[~,data.rf] = sort(data.pos);

for i=1:length(data.strand)
    if(length(unique(data.strand(i,:)))>2)
        error('Several stands at one position.')
    end
end