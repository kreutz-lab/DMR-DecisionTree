function [val, val2] = get(E,prop,arg)
if nargin==1
	disp('-----------------------------')
	disp('Klasse Microarray-Experiment')
	disp(['Name: ',E.name]);
	disp(['data:  ',num2str(size(get(E,'data')))]);
	disp('-----------------------------')
else
	switch lower(prop)

    case {'anzloc'}
        val = sum(~isnan(get(E,'chromloc'))');
    case {'anzna','anznan'}
        val = sum(isnan(get(E,'data')),2);
    case {'antna','antnan'}
        val = sum(isnan(get(E,'data')),2)/size(get(E,'data'),2);
    case {'anznaabsent','anznanabsent'}
        if(isfield(E,'present'))
            val = sum((isnan(get(E,'data')) | ~E.present),2);
        else
            val = ones(get(E,'ngene'),1);
        end
    case {'antnaabsent','antnanabsent'}
        val = get(E,'anznaabsent')./get(E,'na');
    case 'arrays'
        val = E.arrays;
    case 'data'
        ma = get(E,'arrays');
        data = zeros(get(ma{1},'ngene'),get(E,'na'));
        for i=1:get(E,'na')
            data(:,i) = get(ma{i},'data');
        end
		val = data;
    case {'datapresent','dataabsent','datana','datanan'}
        val = get(E,'data');
        present = get(E,'present');
        val(present==0) = NaN;
        case{'genname_short','gen_short','genename_short','genenames_short','gene_short','genes_short'}
            val = get(E,'genenames');
            ind = find(celllength(val)>20);
            for i=1:length(ind)
                val{ind(i)} = val{ind(i)}(1:20);
            end
    case{'genname','gen','genename','genenames','gene','genes'}
        val = get(E,'symbol');
        IDs = get(E,'IDs');
        if(iscell(val))
            nosym = unique([find(celllength(val)==0);strmatch('NA',val,'exact')]);
        else
            nosym = 1:length(IDs);
        end
        
        for i=1:length(nosym)
            val{nosym(i)} = IDs{nosym(i)};
        end
	case 'name'
		val = E.name;
    case {'na','narray','narrays'} % number of arrays
        val = length(E.arrays);
    case {'ngroups','anzg'} % number of groups
        val = length(E.groups);
    case {'ngroupmembers','anzgroupmembers'} % number of groups
        g = E.groups;
        for i=1:length(g)
            lg(i) = length(g{i})
        end
        val = max(lg);
    case {'ngene'} % number of groups
        ma = get(E,'arrays');
        val = get(ma{1},'ngene');
    case {'ngo'} % number of groups
        go = get(E,'go');
        indleer = find(cellempty(go));
        for i=1:length(indleer)
            go{indleer(i)} = 'Data not found';
        end
        val = get(E,'ngene')  - length(find(~cellempty(regexp(go,'Data not found'))));       
    case {'gdata','groupdata'}
        if(~exist('arg','var') | isempty(arg))
            error('get.m: groups number is missing')
        else 
            groupnr = arg;
            data    = get(E,'data');
            groups  = get(E,'groups');
            if(length(arg)>1)
                for i=1:length(groupnr)            
                    val{i} = data(:,groups{groupnr(i)});
                end
            else
                val = data(:,groups{arg});             
            end
        end
    case {'cdata','cfdata'}
        if(~exist('arg','var') | isempty(arg))
            error('get.m: arg is missing')
        else 
            data    = get(E,'data');
            cf      = get(E,'cf',arg);
            lev     = levels(cf);
            for i=1:length(lev)
                val{i}     = data(:,find(cf==lev(i)));
            end
            val2 = lev;
        end
    case {'ddata','dfdata'}
        if(~exist('arg','var') | isempty(arg))
            error('get.m: arg is missing')
        else 
            data    = get(E,'data');
            df      = get(E,'df',arg);
            lev     = levels(df);
            for i=1:length(lev)
                val{i}     = data(:,strcmp(lev{i},df));
            end
            val2 = lev;
        end
        
    case {'hybrname','hybridisierungen','hy'}
        val = E.hybrname;
    case {'groups','gruppen','g'}
        val = E.groups;
    case {'ids','id'}
        val = E.IDs;
        if(size(val,1)==1)
            val = val';
        end

    case {'unigeneid','clusterid','unigene','uid'}
        val = E.unigeneID;
    case 'namen'
        val = emptycells(E.namen);    
        if(size(val,1)==1)
            val = val';
        end
    case 'symbol'
        val = E.symbol;
        if(size(val,1)==1)
            val = val';
        end
    case 'llid'
        val = E.llid;    
    case {'chromosom','chromosome'}
        val = E.chromosome;
    case {'chromloc','chromosomelocation'}
        val = E.chromloc;
    case {'enzym','enzyme'}
        if(isfield(E.annotation,'enzyme'))
            val = E.annotation.enzyme;   
        else
            val = [];
        end
    case {'path','pathway','pathways'}
        if(isfield(E.annotation,'pathways'))
            val = E.annotation.pathways;   
        else
            val = [];
        end
    case {'cn','clonename'}
        val = E.clonename;
    case 'go'
        val = E.go;
    case {'gonames','gonamen'}
        val = E.goNames;
    case 'goids'
        val = E.goIDs;
    case {'df','discretefactors'}
        ma = get(E,'arrays');
        if(~exist('arg') | isempty(arg))
            val = get(ma{1},'df');
        else
            for i=1:length(ma)
                val{i} = get(ma{i},'df',arg);
            end
            val2 = levels(val);
        end
    case {'groupdf','groupdiscretefactors'}
        if(~exist('arg') | isempty(arg))
            error('arg unknown')
        else
            groups = get(E,'groups');
            ma     = get(E,'arrays');
            for i=1:get(E,'ngroups')
                for j=1:length(groups{i})
                    gt{i}{j} = get(ma{groups{i}(j)},'df',arg);
                end
            end
            val = gt;
        end
    case {'cf','continuousfactors'}
        ma = get(E,'arrays');
        if(~exist('arg') | isempty(arg))
            val = get(ma{1},'cf');
        else
            for i=1:length(ma)
                val(i) = get(ma{i},'cf',arg);
            end
            val2 = levels(val);
        end
    case {'groupcf','groupcontinuousfactors'}
        if(~exist('arg') | isempty(arg))
            error('arg unknown')
        else
            groups = get(E,'groups');
            groups = {groups{find(~cellempty(groups))}};
            ma     = get(E,'arrays');
            for i=1:get(E,'ngroups')
                gt{i} = [];
                for j=1:length(groups{i})
                    gt{i} = [gt{i},get(ma{groups{i}(j)},'cf',arg)];
                end
            end
            val = gt;
        end
    case {'meandata','mdata'}
        grouping = arg;
        data = get(E,'data');
        for i=1:length(grouping)
            if(~isempty(grouping{i}))
                val(:,i) = nanmean(data(:,grouping{i}),2);
            else
                val(:,i) = NaN * ones(size(data,1),1);
            end
        end
    case {'mean'}
        val = nanmean(get(E,'data'),2);
    case {'meanpresent'}
        val = nanmean(get(E,'datapresent'),2);
    case {'median'}
        val = nanmedian(get(E,'data')')';
    case {'present'}
        val = E.present;
    case {'sd'}
        val = nanstd(get(E,'data')')';
    case {'sdcall'}
        val = nanstd(get(E,'present')')';
    case {'sdpresent','stdpresent'}
        val = nanstd(get(E,'datapresent')')';
    case {'nan'}
        val = sum(isnan(get(E,'data')'))';
    case {'mm','minmax'}
        val = minmax(get(E,'data')')';
    case {'mmpresent','minmaxpresent'}
        val = minmax(get(E,'datapresent')')';
    case {'mm2','minmax2'}
        val = minmax2(get(E,'data')')';
    case {'genegroups','gengruppen'}
        val = E.genegroups;
    case {'fieldnames'}
        val = fieldnames(E);
    case 'miame'
        if(isfield(E.annotation,'miame'))
            val = E.annotation.miame;   
        else
            val = [];
        end
    case {'rfnadel','reihenfolgepronadel'}
        val = GetReihenfolgeProNadel(E);
	otherwise
        try
            val = E.(prop);
        catch
            try
                val = E.container.(prop);
            catch
    		    error('maExperiment/get.m: Property ',prop,' unknown.');
            end
        end
	end
end

if(~exist('val2','var'))
    val2 = [];
end
