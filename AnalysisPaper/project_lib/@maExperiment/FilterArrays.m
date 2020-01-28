% Extrahiert ein Subset von Arrays


function Efilter = FilterArrays(E,filter)
hybrname = get(E,'hybrname');
Efilter = E;
na      = get(E,'na');
ma      = get(E,'arrays');


if(iscell(filter))
    drin = [];
    for i=1:length(filter)
        tmp = strmatch(filter{i},hybrname,'exact');
        if(~isempty(tmp))
            drin = [drin,tmp];
        end
    end
    filter = drin;
end

if( sum(filter==0  | filter==1) == length(filter)  & sum(size(filter)>2))
    if(length(filter) == get(E,'na'))
        filter = find(filter);  % wandelt boolsche Varialble in Indizes um.
    else
        error('FilterGene.m: length(filter) must be equal number of arrays.')
    end
end


anz = 0;

maneu = cell(0);
hybrneu = cell(0);
for i=1:length(filter)
    anz = anz+1;
    maneu{anz} = ma{filter(i)};
    hybrneu{anz} = hybrname{filter(i)};
end

Efilter.arrays   = maneu;
Efilter.hybrname = hybrneu;

if(~isempty(E.present))
    Efilter.present = E.present(:,filter);
end

gr = get(E,'groups');
for i=1:length(gr)
    gr{i} = intersect(gr{i},filter);
    for j=1:length(gr{i})
        gr{i}(j) = gr{i}(j) - (gr{i}(j) - sum(filter<=gr{i}(j)));
    end
end
Efilter.groups = gr;

% try
%     Efilter.groups = FindGroups(hybrneu, length(get(E,'groups')) );
% catch
%     disp('FilterArrays.m: Gruppierung fehlgeschlagen')
% end

cont = get(E,'container');
% fn = fieldnames(cont);
% for i=1:length(fn)
%     if(size(cont.(fn{i}),2)==na)
%         cont.(fn{i}) = cont.(fn{i})(:,filter);
%     end
% end
cont = FilterStruct(cont,filter,na,2);

% if(isfield(cont,'pdata'))
%     cont.pdata.SampleNr
%     cont.pdata = FilterStruct(cont.pdata,filter,na,2);
%     cont.pdata.SampleNr
% end
Efilter = set(Efilter,'container',cont);

