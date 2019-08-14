function dats = wgbss_CalculateDataFeatures(dats)
% if ~exist('thresh','var') || isempty(thresh)
%     calcthresh = true
% else
%     calcthresh = false;
% end

for i=1:length(dats)    
%     if max(dats(i).freqT)>1
%         dats(i).freqC = 1-0.01*dats(i).freqT;
%     end
%     if calcthresh
        threshObjFun = quantile(dats(i).smooth.freqC(:),[.2,.8]);
%         threshObjFun(1) = max(threshObjFun(1),1-threshObjFun(2));
%         threshObjFun(2) = min(threshObjFun(2),1-threshObjFun(1));

%         thresh = mean(quantile(dats(i).smooth.freqC,[.2,.8]));
        dats(i).thresh_ismeth = threshObjFun;
%     end
    

    dats(i).isMeth = logical(size(dats(i).pos));
    indicator = FindCloser(dats(i).smooth.freqC(dats(i).rf,1)<=threshObjFun(1),  dats(i).smooth.freqC(dats(i).rf,1) >=threshObjFun(2),  dats(i).pos(dats(i).rf))';
    dats(i).isMeth(dats(i).rf) =  indicator>= 1.5;
    
    if sum(dats(i).isMeth)>2    
        dats(i).Meth.mean = nanmean(dats(i).freqC(dats(i).isMeth));
        dats(i).Meth.sd = nanstd(dats(i).freqC(dats(i).isMeth));
        dats(i).Meth.meanNread = nanmean(dats(i).nread(dats(i).isMeth));
        dats(i).Meth.sdNread = nanstd(dats(i).nread(dats(i).isMeth));
    else
        dats(i).Meth.mean = 0;
        dats(i).Meth.sd = 0;     
        dats(i).Meth.meanNread = 0;
    end
    
    if sum(~dats(i).isMeth)>2    
        dats(i).NonMeth.mean = nanmean(dats(i).freqC(~dats(i).isMeth));
        dats(i).NonMeth.sd = nanstd(dats(i).freqC(~dats(i).isMeth));    
        dats(i).NonMeth.meanNread = nanmean(dats(i).nread(~dats(i).isMeth));
        dats(i).NonMeth.sdNread = nanstd(dats(i).nread(~dats(i).isMeth));
    else
        dats(i).NonMeth.mean = 0;
        dats(i).NonMeth.sd = 0;     
        dats(i).NonMeth.meanNread = 0;
    end
    
    pos = dats(i).pos(dats(i).rf);
    
    % Abstand CpGs
    diffpos = diff(pos);
    dats(i).meanLog10DiffPos = nanmean(log10(diffpos(diffpos>0)));
    dats(i).sdLog10DiffPos = nanstd(log10(diffpos(diffpos>0)));    

    indup = find(diff(dats(i).isMeth(dats(i).rf))==1);
    inddown = find(diff(dats(i).isMeth(dats(i).rf))==-1);
    indboth = union(indup,inddown);

    if isempty(indup) || isempty(inddown)
        save('error_noup_or_nodown')
    end
    
%     nend = floor(length(indboth)/2)*2;  % gerade Zahl
%     uend = nend-1;  % ungerade Zahl
    if ~isempty(indup) && ~isempty(inddown) && indup(1)>inddown(1) % erst hoch
        % Länge Methylierung
        dats(i).meanLog10MethLengthPos = nanmean(log10( pos(indboth(2:2:end))-pos(indboth(1:2:(end-1)))  ));
        dats(i).sdLog10MethLengthPos   =  nanstd(log10( pos(indboth(2:2:end))-pos(indboth(1:2:(end-1))) ));
        
        % Länge Demethylierung
        dats(i).meanLog10NonMethLengthPos = nanmean(log10( pos(indboth(3:2:end))-pos(indboth(2:2:(end-1)))  ));
        dats(i).sdLog10NonMethLengthPos = nanstd(log10(    pos(indboth(3:2:end))-pos(indboth(2:2:(end-1)))  ));
    else % erst runter
        % Länge Methylierung
        dats(i).meanLog10MethLengthPos = nanmean(log10(  pos(indboth(3:2:end))-pos(indboth(2:2:(end-1)))  ));
        dats(i).sdLog10MethLengthPos   = nanstd(log10(   pos(indboth(3:2:end))-pos(indboth(2:2:(end-1)))  ));
        
        % Länge Demethylierung
        dats(i).meanLog10NonMethLengthPos = nanmean(log10(  pos(indboth(2:2:end))-pos(indboth(1:2:(end-1)))  ));
        dats(i).sdLog10NonMethLengthPos   = nanstd(log10(  pos(indboth(2:2:end))-pos(indboth(1:2:(end-1)))  ));
    end
end


