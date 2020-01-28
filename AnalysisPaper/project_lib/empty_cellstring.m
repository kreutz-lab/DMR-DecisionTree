% cs = empty_cellstring(dim1,dim2,value)
% 
% 1. Zeile code, die eine Zelle von leeren String '' erzeugt.
%   cs = strrep(num2strArray(ones(dim1,dim2)),'1','');

function cs = empty_cellstring(dim1,dim2,value)
if(nargin==2)
    value = '';
end
% cs = strrep(num2strArray(ones(dim1,dim2)),'1',value);
cs = rep({value},[dim1,dim2]);
