% tabs = tablesMatchOrder(tabs, orderfield)
% 
%   tabs    cell of tables
% 
%   This function produces tables of the same row order and same size. If
%   required, new rows are added with NaN or 'NA' entries.

function tabs = tablesMatchOrder(tabs, orderfield)

nrow = NaN(1,length(tabs));
for t=1:length(tabs)
    if length(unique(tabs{t}.(orderfield)))< length(tabs{t}.(orderfield))
        error('tabs{t}.%s is have non-unique entries',t,orderfield)
    end
    nrow(t) = size(tabs{t},1);
end
[nrowmax,iref] = max(nrow);
iref = iref(1); % the order in this subset will be used as reference

% create a one-row table with NaNs
emptyline = table2cell(tabs{1}(1,:));
emptyline(1,cellfun(@isnumeric,emptyline)) = {NaN};
emptyline(1,cellfun(@ischar,emptyline)) = {'NA'};
emptyline = cell2table(emptyline,'VariableNames',tabs{1}.Properties.VariableNames);

for j=setdiff(1:length(tabs),iref)
    
    [~,ia,ib] = intersect(tabs{iref}.(orderfield),tabs{j}.(orderfield),'stable');
    tmp = rep(emptyline,[nrowmax,1]);
    tmp(ia,:) = tabs{j}(ib,:);
    tmp.Properties.Description = tabs{j}.Properties.Description;
    tabs{j} = tmp;
end
