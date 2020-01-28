function [rang,F1,rangSorted] = ImageBest(names_meth,names_context,ass2x,ass2y,xval)
% names_context = lower(names_context);
    addpath('Rlink')

methlev = levels(names_meth);
names_context = strrep(names_context,'chg','CHG');
contextlev = levels(names_context);

if length(xval)>1
    for ix=1:length(xval)
        [rang(:,:,ix),ass2(:,:,ix)] = ImageBest(names_meth,names_context,ass2x,ass2y,xval(ix));
    end
    F1 = mean(ass2,3);
    rang = mean(rang,3);

    E = maExperiment(rang);
    erg = ClusterLarge(E);
    system('rm Cluster*'); % Heatmap etc not required
    
    erg.orderChips(erg.orderChips==strmatch('MethylScore',methlev,'exact'))=[]; % unpublished and no compliance from authors

    rangSorted = rang(erg.order,erg.orderChips);
    imagesc_nan(rangSorted,colormap('cool'));
    title('Rank');
    colorbar
    xticklabel_rotate(1:size(rangSorted,2),90,methlev(erg.orderChips))
    set(gca,'YTick',1:size(rang,1),'YTickLabel',ReplaceSampleNames(contextlev(erg.order)),'FontSize',6);
    saveas(gcf,'ImageBest');
    print -dpng ImageBest

    figure
    imagesc_nan(F1(erg.order,erg.orderChips),colormap('cool'));
    title('F1-score');
    colorbar
    xticklabel_rotate(1:size(rangSorted,2),90,methlev(erg.orderChips))
    set(gca,'YTick',1:size(rangSorted,1),'YTickLabel',ReplaceSampleNames(contextlev(erg.order)),'FontSize',6);
    saveas(gcf,'ImageF1');
    print -dpng ImageF1

    figure
    set(gca,'FontSize',14)
    boxplot(rangSorted,'plotstyle','compact','Labels',methlev(erg.orderChips))
    ylabel('Rank')
%     publ(2)
%     xticklabel_rotate(1:size(rang,2),90,methlev)
    set(gca,'FontSize',14)
    saveas(gcf,'BoxplotBest');
    print -dpng BoxplotBest

    figure
    set(gca,'FontSize',14)
    boxplot(F1(erg.order,erg.orderChips),'plotstyle','compact','Labels',methlev(erg.orderChips))
    ylabel('F1-score')
    set(gca,'FontSize',14)
%     publ(2)
%     xticklabel_rotate(1:size(rang,2),90,methlev)
    saveas(gcf,'BoxplotF1');
    print -dpng BoxplotF1
else
    
    
    rang = NaN(length(contextlev),length(methlev));
    F1 = NaN(length(contextlev),length(methlev));
    
    for c=1:length(contextlev)
        ass2 = NaN(size(methlev));
        for m=1:length(methlev)
            imeth = strmatch(methlev{m},names_meth,'exact');
            icon = strmatch(contextlev{c},names_context,'exact');
            
            irow = intersect(imeth,icon);
            if length(irow)>1
                irow
                error('multiple match')
            elseif isempty(irow)
                warning('No match for %s %s\n',methlev{m},contextlev{c})
            else
                fprintf('%s & %s found\n',methlev{m},contextlev{c})
            end
            ix = find(ass2x(irow,:)<=xval);
            
            if ~isempty(ix)
                ass2(m) = ass2y(irow,ix(end));
            end
        end
        %     [~,rf] = sort(-ass2);
        % ass2
        F1(c,:) = ass2;
        r = rankasgn_fast(ass2);
        r(isnan(ass2)) = NaN;
        r = r-min(r); % min =0
        if max(r)>0
            r = r./max(r); % between 0 and 1
            rang(c,:) = r;
        end
    end
end


