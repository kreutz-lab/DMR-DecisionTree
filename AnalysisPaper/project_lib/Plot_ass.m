function Plot_ass(ass,file,tit,X,xnames)
if ~exist('tit','var')
    tit = '';
end
if ~exist('X','var') || isempty(X)
    X = ones(length(ass),1);
end
if ~exist('xnames','var') || isempty(xnames)
    xnames = cell(0);
end

fn = fieldnames(ass(1));
close all
figure
set(gcf,'Position',1e3*[0.0081    0.2050    1.7004    0.6529])
subx = ceil(sqrt(length(fn)));
suby = ceil(length(fn)/subx);
for i=1:length(fn)
    subplot(suby,subx,i)
    momass = [ass.(fn{i})];
    for ix=1:size(X,2)
        tmp = momass(X(:,ix)==1);
        if size(tmp,1) ==1
            tmp = tmp';
        end
        plotdat(:,ix) = tmp;
    end
    if isempty(X)
        hist(plotdat,50)
    else 
        hist(plotdat,10)
    end
    
    title(fn{i})
    set(gca,'FontSize',7)
end
if ~isempty(tit)
    suptitle(tit)
end
if ~isempty(xnames)
    legend(str2label(xnames))
end
PrintToPng(gcf,file)
