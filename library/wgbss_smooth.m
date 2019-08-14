function [ys2,thresh] = wgbss_smooth(y,x,doplot)

if ~exist('doplot','var') || isempty(doplot)
    doplot = 0;
end
global fastSmooth

if fastSmooth==1
    spans = 10;
else
    spans = [10,20,50];
end

if size(y,1)>1 && size(y,2)>1
    ys2 = NaN(size(y));
    thresh = NaN(size(y,2),2);
    for i=1:size(y,2)  
        [ys2(:,i),thresh(i,:)] = wgbss_smooth(y(:,i),x,doplot);
    end    
else

    if size(y,1)>1
        dotranspose = 1;
        y = y'; % make row
    else 
        dotranspose = 0;
    end
    ys = NaN(length(spans),length(y));
    
    for i=1:length(spans)
        if fastSmooth
            ys(i,:) = smooth(y,spans(i),'moving');
        else
            ys(i,:) = smooth(y,spans(i),'lowess');
        end        
    end
    
    thresh = quantile(ys(:),[.2,.8]);
    if max(ys(:))<=1  % fractions in range [0,1]
        thresh(2) = min(thresh(2),.8);
        thresh(1) = max(thresh(1),.2);
    elseif abs(max(ys(:))-100) <0.1  % percent in range [0,100]
        thresh(2) = min(thresh(2),80);
        thresh(1) = max(thresh(1),20);
    end
        
%     bol(1,:) = min(ys,[],1)<thresh(1);
    bol(1,:) = y<=thresh(1);
    bol(3,:) = y>=thresh(2);
%     bol(3,:) = max(ys,[],1)>thresh(2) & ~bol(1,:) ;
    bol(2,:) = median(ys,1)>thresh(1) & median(ys,1)<thresh(2) & ~bol(1,:) & ~bol(3,:);
        
    ys2 = NaN(1,size(ys,2));
    
    ys2(bol(1,:)) = min(ys(:,bol(1,:)),[],1);
    ys2(bol(3,:)) = max(ys(:,bol(3,:)),[],1);
        
    indicator = FindCloser(bol(1,:),bol(3,:),x);
    
    bol1_closer = bol(2,:) & indicator ==1;
    bol3_closer = bol(2,:) & indicator ==2;
    
    ys2(bol1_closer) = min(ys(:,bol1_closer),[],1);
    ys2(bol3_closer) = max(ys(:,bol3_closer),[],1);
        
    % ys2(bol(2,:)) = median(ys(:,bol(2,:)),1);
    ys2(bol(3,:)) = max(ys(:,bol(3,:)),[],1);
    ys2(bol(1,:)) = min(ys(:,bol(1,:)),[],1);
    
    bol_nan = isnan(ys2);
    if sum(bol_nan)>0
        ys2(bol_nan)   = mean(ys(:,bol_nan),1);
    end    
    
    if dotranspose
        ys2 = ys2';
    end
    
    if doplot~=0
        h2=plot(x,y,'k.');
        hold on
        if doplot==1
            h1=plot(x,ys');
            h3=plot(x,ys2,'m','LineWidth',2);
        elseif doplot==3
            h3=plot(x,ys2,'m','LineWidth',1);
        end
        
        threshObjFun = quantile(ys2,[.2,.8]);
%         threshObjFun(1) = max(threshObjFun(1),1-threshObjFun(2));
%         threshObjFun(2) = min(threshObjFun(2),1-threshObjFun(1));
    threshObjFun(2) = min(threshObjFun(2),mean([.8,threshObjFun(2)]));
    threshObjFun(1) = max(threshObjFun(1),mean([.1,threshObjFun(1)]));

        indnon = find(ys2<=threshObjFun(1));
        indmeth = find(ys2>=threshObjFun(2));  
        indboth = union(indnon,indmeth);
        
        ymeth = NaN(size(x));
        ymeth(indmeth) = 1;
        ynon = NaN(size(x));
        ynon(indnon) = 0;
        if doplot==1
            h4 = plot(x,ymeth+.02, '-', 'Color',.7*ones(1,3),'LineWidth',2);
            h4 = plot(x,ynon-.02, '-', 'Color',.7*ones(1,3),'LineWidth',2);
        end
%         ymeth2 = ones(size(indboth));
%         [~,ia] = intersect(indboth,indnon);
%         ymeth2(ia) = NaN;
%         h5=plot(x(indboth),ymeth2+.05,'r','LineWidth',3);
%         
%         ynon2 = zeros(size(indboth));
%         [~,ia] = intersect(indboth,indmeth);
%         ynon2(ia) = NaN;
%         h5 = plot(x(indboth),ynon2-.05,'r','LineWidth',3);

        ismeth = round(FindCloser(ys2<=threshObjFun(1),ys2>=threshObjFun(2),x)-1);        
        
        if doplot==1 || doplot==2
            h5 = plot(x,ismeth*1.1-0.05,'Color','r','LineWidth',1);
        elseif doplot==3
            h5 = plot(x,ismeth*1.1-0.05,'Color','r','LineWidth',1);
        end
        
        set(gca,'FontSize',14)
        xlabel('pos')
        ylabel('freqC')
        tmp = strcat('span=',num2strArray(spans));
        if doplot==1
            set(legend([h2;h1;h3;h4;h5],'Data',tmp{:},'Used smoothing','Individual calls','Call regions'),'FontSize',12)
        end
        axis tight
%         xlim([525240    664360])
%         ylim([-.1,1.1])

    end
end

global doSaveSmoothing
if doSaveSmoothing
    save wgbss_smooth
end
