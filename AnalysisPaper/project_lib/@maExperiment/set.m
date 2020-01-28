% set(M,proptery, value, argument)
% argument ist z.B. das Feld bei M.continuousFactors

function E = set(E,varargin)
property_argin=varargin;
while length(property_argin) >=2
    prop = property_argin{1};
    val  = property_argin{2};
    if(length(property_argin)>2)
        arg  = property_argin{3};
    else 
        arg = [];
    end
    property_argin=property_argin(4:end);
    
    switch lower(prop)
    case 'name'
        E.name = val;
    case 'arrays'
        E.arrays = val;
    case {'chromosom','chromosome'}
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.chromosome = val;    
        else
            error('Length incorrect.')
        end
    case {'chromloc','chromosomelocation'}
        if(size(val,1)==size(get(E,'data'),1) | isempty(val))
            E.chromloc = val;    
        else
            error('Length incorrect.')
        end
    case {'cn','clonename'}
        if(size(val,1)==size(get(E,'data'),1) | isempty(val))
            E.clonename = val;    
        else
            error('Length incorrect.')
        end
    case {'cf','continuousfactors'}
        if(isempty(val))
            ma = get(E,'arrays');
            for i=1:get(E,'na')            
                ma{i} = set(ma{i},'cf',[]);
            end
            E = set(E,'arrays',ma);
        elseif(length(val)~= get(E,'na'))
            error('number of time points is invalid')
        elseif(~exist('arg') | isempty(arg))
                error('arg unknown')
        else
            ma = get(E,'arrays');
            for i=1:get(E,'na')            
                ma{i} = set(ma{i},'cf',val(i), arg);
            end
            E = set(E,'arrays',ma);
        end
    case 'data'
        for i=1:size(val,2)
            E.arrays{i} = set(E.arrays{i},'data',val(:,i));
        end
    case {'enzym','enzyme'}
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.annotation.enzyme = val;    
        else
            error('Length incorrect.')
        end
    case 'go'
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.go = val;    
        else
            error('Length incorrect.')
        end
    case 'gonames'
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.goNames = val;    
        else
            error('Length incorrect.')
        end
    case 'goids'
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.goIDs = val;    
        else
            error('Length incorrect.')
        end
    case {'hybrname','hybridisierungen','hy'}
        E.hybrname = val;
        for i=1:length(E.arrays)
            E.arrays{i} = set(E.arrays{i},'name',E.hybrname{i});            
        end
    case {'ids','id'}
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.IDs = val;
        else
            error('Length incorrect.')
        end
    case 'llid'
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.llid = val;    
        else
            error('Length incorrect.')
        end
    case 'namen'
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.namen = val;    
        else
            error('Length incorrect.')
        end
    case {'path','pathway','pathways'}
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.annotation.pathways = val;
        else
            error('Length incorrect.')
        end
    case 'symbol'
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.symbol = val;    
        else
            error('Length incorrect.')
        end
    case {'unigeneid','clusterid','unigene','uid'}
        if(length(val)==size(get(E,'data'),1) | isempty(val))
            E.unigeneID = val;
        else
            error('length incorrect');
        end
    case {'groupcf'} % z.B. die Zeiten für die Chips einer Gruppe 
        if(length(val)~= get(E,'ngroups'))
            error('invalid number of grouptimes')
        elseif(~exist('arg') | isempty(arg))
            error('arg unknown')
        else
            ts = [];
            groups = get(E,'groups');
		    for i=1:length(groups)
                if(~isempty(groups{i}))
                    ts(groups{i}) = val(i);
                end
		    end
    		E = set(E,'cf',ts, arg);
        end

    case {'df','treatments'}
        if(isempty(val))
            ma = get(E,'arrays');
            for i=1:get(E,'na')            
                ma{i} = set(ma{i},'df',[]);
            end
            E = set(E,'arrays',ma);
        elseif(length(val)~= get(E,'na'))
            error('number of treatments is invalid')
        elseif(~exist('arg') | isempty(arg))
                error('arg unknown')
        else
            ma = get(E,'arrays');
            if(isnumeric(val))
                val = array2stringcell(val);
            end
            for i=1:get(E,'na')            
                ma{i} = set(ma{i},'df',val{i},arg);
            end
            E = set(E,'arrays',ma);
        end
        
    case {'groupdf','groupdiscretefactors'} % z.B. die Behandlung der Probe der Chips einer Gruppe 
        if(length(val)~= get(E,'ngroups'))
            error('invalid number of grouptreatments')
        elseif(~exist('arg') | isempty(arg))
            error('arg unknown')
        else
            ts = cell(0);
            groups = get(E,'groups');
    		for i=1:length(groups)
                if(~isempty(groups{i}))
                    for j=1:length(groups{i})
                        ts{groups{i}(j)} = val{i};
                    end
                end
            end
    		E = set(E,'df',ts, arg);
        end
   case {'genegroups','gengruppen'}
        E.genegroups = val;
   case {'addgenegroups','addgg'}
        tmp = E.genegroups;
        tmp{end+1} = val;
        E.genegroups = tmp;
   case {'groups','gruppen','group'}
        E.groups = val; 
   case {'miame'}
        E.annotation.miame = val; 
   otherwise
       try
           E.(prop) = val;
       catch
           E.container.(prop) = val;
           disp(['maExperiment/set.m: Property ',prop,' unknown, in Container gesteckt.']);
       end
   end
end
