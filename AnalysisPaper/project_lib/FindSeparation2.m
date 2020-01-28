% res = FindSeparation(targetfield, tab1, tab2, rankmeth, ignorevars)
%
% Example:
%   res = FindSeparation('F1', tabs([1,7]), 1, {'method'}, 'context')


function res = FindSeparation2(tabs, targetfield, ignorevars, orderfield)
if ~exist('ignorevars','var') || isempty(ignorevars)
    cell(0);
elseif ischar(ignorevars)
    ignorevars = {ignorevars};
end

tabs = tablesMatchOrder(tabs, orderfield);

target1 = tabs{1}.(targetfield);
target2 = tabs{2}.(targetfield);

% remove
[~,drin] = setdiff(tabs{1}.Properties.VariableNames,[ignorevars,{targetfield}],'stable');
for i=1:length(tabs)
    tabs{i} = tabs{i}(:,drin);
end

if length(tabs)~=2
    error('Current implementation is only for two tables')
end

res = struct;
for j=1:size(tabs{1},2)
    best1 = find(target1>target2);
    best2 = find(target1<target2);

    x = table2array(tabs{1}(best1,j));
    y = table2array(tabs{1}(best2,j));    
    x2 = table2array(tabs{2}(best1,j));
    y2 = table2array(tabs{2}(best2,j));    
    
    if sum(abs(x-x2))>0
        tabs{1}.Properties.VariableNames{j}
        x-x2
        error('tab and tab2 should have same predictor values (x)')
    end
    if sum(abs(y-y2))>0
        tabs{1}.Properties.VariableNames{j}
        y-y2
        error('tab and tab2 should have same predictor values (y)')
    end
    
    if ~isempty(x) && ~isempty(y)
        res.p_ranksum(1,j) = ranksum(x,y);
        res.p_optimal(1,j) = ranksum(ones(size(x)),zeros(size(y)));
    else
        res.p_ranksum(1,j) = 1;
        res.p_optimal(1,j) = 1;
    end
    
    res.name{j} = tabs{1}.Properties.VariableNames{j};
end

res.isoptimal = res.p_ranksum<=(res.p_optimal+1e-4);


