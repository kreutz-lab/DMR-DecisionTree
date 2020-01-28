% tabs = table2subset(tab,subfield)
% tabs = table2subset(tab,subfield,orderfield)
% 
% This function splits a table according to the levels of tab.(subfield)
% 
%   orderfield  if specified, all output tables subsets have the same
%               number of rows and the same order. Ordering is done by the
%               values in  tab.(orderfield). These values have to be unique
%               in each subset
% 
% 
% Examples:
% tabs = table2subset(simDatF,'method')
% tabs = table2subset(simDatF,'method','context')

function tabs = table2subset(tab,subfield,orderfield)
if ~exist('orderfield','var')
    orderfield = [];
end

isub = strmatch(subfield,tab.Properties.VariableNames,'exact');
if length(isub)~=1
    subfield
    tab.Properties.VariableNames
    error('subfield does not match uniquely.')
end
subvals = tab.(subfield);
levs = levels(subvals);

tabs = cell(size(levs));
for i=1:length(levs)
    ind = ismember(tab.(subfield),levs(i));    
    tabs{i} = tab(find(ind),:);
    if iscell(levs(i))
        sprintf('%s = %s',subfield,levs{i})
        tabs{i}.Properties.Description = sprintf('%s = %s',subfield,levs{i});
    end
end

if ~isempty(orderfield) % enforce the same order according to tab.(orderfield)
    tabs = tablesMatchOrder(tabs, orderfield);
end
