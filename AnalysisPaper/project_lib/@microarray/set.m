% set(M,proptery, value, argument)
% argument ist z.B. das Feld bei M.continuousFactors

function M = set(M,varargin)

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
        M.name = val;
    case 'data'
%         if(length(val) == length(M.data))
            M.data = val;
%         else
%             error(['microarray/set.m: length(val)=',num2str(length(val)),' ~= length(M.data)=',num2str(length(M.data))])
%         end
        
    case {'df','discretefactors'}
        if(~exist('arg'))
            error('arg unknown')
        elseif(isempty(arg))
            M.discreteFactors = [];
        else
            M.discreteFactors = setfield(M.discreteFactors,arg, val);
        end
    case {'cf','continuousfactors'}
        if(~exist('arg'))
            error('arg unknown')
        elseif(isempty(arg))
            M.continuousFactors = [];
        else
            M.continuousFactors = setfield(M.continuousFactors,arg, val);
        end
        
    otherwise
        error('microarray/set.m: Property unknown.');
    end
end
