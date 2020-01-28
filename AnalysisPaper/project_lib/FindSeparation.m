% res = FindSeparation(targetfield, tab1, tab2, rankmeth, ignorevars)
%
% Example:
%   res = FindSeparation('F1', tabs([1,7]), 1, {'method'}, 'context')


function res = FindSeparation(targetfield, tabs, rankmeth, ignorevars, orderfield)
if ~exist('rankmeth','var') || isempty(rankmeth)
    rankmeth = 1;
end
if ~exist('ignorevars','var') || isempty(ignorevars)
    cell(0);
elseif ischar(ignorevars)
    ignorevars = {ignorevars};
end

tabs = tablesMatchOrder(tabs, orderfield);
scores = NaN(size(tabs{1},1),length(tabs));
for i=1:length(tabs)
    scores(:,i) = tabs{i}.(targetfield);
end

if rankmeth<0
    rankmeth = -rankmeth;
    scores = -scores;
end

ranks = NaN(size(scores));
for i=1:size(scores,1)
    switch rankmeth
        case 1 % default
            ranks(i,:) = Rankasgn(scores(i,:));
        case 2
            notnan = find(~isnan(scores(i,:)));
            r = Rankasgn(levs(i));
            r = (r-1)/(length(notnan)-1)*2-1; % scaling to -1,1 for ~isnan
            ranks(i,notnan) = r(notnan);
            
        otherwise
            rankmeth
            error('rankmeth unknown.');
    end
end

% remove
[~,drin] = setdiff(tabs{1}.Properties.VariableNames,[ignorevars,{targetfield}],'stable');
for i=1:length(tabs)
    tabs{i} = tabs{i}(:,drin);
end

res = struct;

if length(tabs)~=2
    error('Current implementation is only for two tables')
end

% res.p_ranksum = NaN(length(rs),size(preds,2));
[~,isbest] = min(ranks,[],2);

for j=1:size(tabs{1},2)
%     try
    x = table2array(tabs{1}(isbest==1,j));
    y = table2array(tabs{1}(isbest==2,j));
    res.p_ranksum(1,j) = ranksum(x,y);
    res.p_optimal(1,j) = ranksum(ones(size(x)),zeros(size(y)));
    res.name{j} = tabs{1}.Properties.VariableNames{j};
%     catch
%         isbest
%     end
end

%     res.isoptimal(i,:) = (res.p_ranksum(i,:)==ranksum(zeros(1,sum(isbest)),ones(1,sum(~isbest)))); % perfect separation, same sample size



