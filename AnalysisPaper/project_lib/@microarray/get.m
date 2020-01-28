function val = get(M,prop,arg);
if nargin==1
	disp('-----------------------------')
	disp('Klasse Microarray')
	disp(['Name: ',get(M,'name')]);
	disp(['data:  ',num2str(size(get(M,'data')))]);
    if(isstruct(M.discreteFactors))
        f = fieldnames(M.discreteFactors);
        nfield = length(f);
        disp('discrete   Factors:')        
        for i=1:nfield
            disp(f{i});
        end
    end        
    if(isstruct(M.continuousFactors) )
        f = fieldnames(M.continuousFactors);
        nfield = length(f);
        disp('continuous Factors:')        
        for i=1:nfield
            disp(f{i});
        end
    end    
	disp('-----------------------------')
else
	switch prop
	case 'name'
		val = M.name;
	case 'data'
		val = M.data;
	case {'ngene','ng'}
		val = size(M.data,1);
    case {'df','discreteFactors'}
        if(~exist('arg') | isempty(arg))
            if(isstruct(M.discreteFactors))
                val = fieldnames(M.discreteFactors);
            else
                val = [];
            end
        else
            val = getfield(M.discreteFactors,arg);
        end
    case {'cf','continuousFactors'}
        if(~exist('arg') | isempty(arg))
            if(isstruct(M.continuousFactors))
                val = fieldnames(M.continuousFactors);
            else
                val = [];
            end
        else
            val = getfield(M.continuousFactors,arg);
        end
	otherwise
		error('microarray/get.m: Property unknown.');
	end
end
