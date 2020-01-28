% tree = fitTreeToEachMethod(X,y,xnames,varargin)

function [tree,handles] = fitTreeToEachMethod(X,y,xnames,methnames, prefix, varargin)

tree_args = {'MinParentSize',2,...
    'MaxNumSplits',1,...
    'Surrogate',3,...
    'MinLeafSize',1,...
    'HyperparameterOptimizationOptions',struct('SaveIntermediateResults',true,'Verbose',2,'ShowPlots',true),...
    };

for v= 1:2:length(varargin)
    tree_args{end+1} = varargin{v};
    tree_args{end+1} = varargin{v+1};
end    

% if method is a predictor
imeth = strmatch('method',lower(xnames));
if length(imeth)>1
    error('length(imeth)>1')
elseif length(imeth)==1   
    xmeth = X(:,imeth);    
    levm = levels(xmeth);
    X(:,imeth) = [];
    xnames(imeth) = [];
else
    error('method should be one column in X labelled as xnames=''method''.');
end



tree = cell(size(levm));
handles = cell(0);
for i=1:length(levm)
    ia = ismember(xmeth,levm(i));
    
    tree{i} = fitrtree(X(ia,:),y(ia),tree_args{:},'PredictorNames',xnames);
    [handles{i},axhandle] = viewtree(tree{i});
    
    htit = get(axhandle,'Title');
    titstr = methnames{i};
    set(htit,'String',titstr,'Position',[0.5,.95,0]);
    
    warning('off','stats:classreg:learning:regr:CompactRegressionTree:get:Risk');
    s = struct(tree{i});
    if s.NumNodes==3
        
        labels = plotTree(handles{i});
        
        spl = strsplit(strtrim(labels{1}),' ');
        ifn = strmatch(spl{1},xnames,'exact');
        evstr = ['sum(X(ia,',num2str(ifn),')',spl{2},spl{3},');'];
        nstr = ['n=',num2str(eval(evstr))];
        text(-1-.2,-.4,nstr,'HorizontalAlignment','left','FontSize',19);
        
        spl = strsplit(strtrim(labels{2}),' ');
        ifn = strmatch(spl{1},xnames,'exact');
        evstr = ['sum(X(ia,',num2str(ifn),')',spl{2},spl{3},');'];
        nstr = ['n=',num2str(eval(evstr))];
        text(1+.2,-.4,nstr,'HorizontalAlignment','right','FontSize',19);
        
        
        title(titstr,'FontSize',20)
        print(sprintf('%sTree_%i.png',prefix,i),'-dpng','-r600');
        
        close(handles{i});
    else     
        set(handles{i},'Position',[0  255  1200  315])
        print(handles{i},sprintf('%sTree_%i.png',prefix,i),'-dpng','-r600');
        close(handles{i});
    end
end

