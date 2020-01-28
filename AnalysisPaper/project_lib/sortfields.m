% sout = sortfields(sin)
% 
%   Sortiert die Felder eines Structs alphabetisch.
% 
%   sin     struct
% 
%   sout    struct

function sout = sortfields(sin)
fn = sort(fieldnames(sin));

for i=1:length(fn)
    sout.(fn{i}) = sin.(fn{i});
end
