% s = FilterStruct(s,filter,ndat,dim)
% 
%   Alle Felder im Struct werden gefiltert, wenn die Längen entlang der 1.
%   oder 2. Dimension identisch ndat ist.
% 
%   filter      Indizes, die verbleiben
%   ndat        Länge der Matrizen, Array, Zellen (entlang der zu filternden
%               Dimension)
%   dim         Falls dim existiert, so wird nur entlang dieser Dimension
%               gefiltert. 
%               Bisher nur dim=[], dim=1, dim=2 zulässig.

function s = FilterStruct(s,filter,ndat,dim)
if(nargin==0)
    disp('s = FilterStruct(s,filter,ndat,dim)')
    return
end
if(~exist('dim','var') | isempty(dim))
    dim = [];
elseif(dim~=1 & dim~=2)
    error('FilterStruct.m: Bisher nur dim=[], dim=1, dim=2 zulässig. ')
end
if(~exist('ndat','var') | isempty(ndat))
    fn = fieldnames(s);
    ndat = length(s.(fn{1}));
end

if(isstruct(s))
    fn = fieldnames(s);
    for i=1:length(fn)
        fval = s.(fn{i});
        if(isstruct(fval))
            fval = FilterStruct(fval,filter,ndat,dim);
            s.(fn{i}) = fval;
        elseif(size(fval,1)==ndat & (isempty(dim) | dim==1) )       % Filter über die 1. Dimension der Matrix/Zelle
            s.(fn{i})=Filter2dObjekt(getfield(s,fn{i}),filter);
        elseif(size(fval,2)==ndat & (isempty(dim) | dim==2))        % Filter über die 2. Dimension der Matrix/Zelle
            s.(fn{i}) = Filter2dObjekt(s.(fn{i})',filter)';
        else
%             disp(['FilterStruct.m: Cannot filter Property ',fn{i}]);
        end
    end
end
