function dats = wgbss_SmoothDat(dats,doplot)

if ~exist('doplot','var') || isempty(doplot)
    doplot = 0;
end
global doSaveSmoothing


for i=1:length(dats)
    if ~isfield(dats(i),'freqC')
        if max(dats(i).freqT)>1
            dats(i).freqC = 100 - dats(i).freqT;
        else
            dats(i).freqC = 1 - dats(i).freqT;
        end
    end
    
    [dats(i).smooth.freqC(dats(i).rf),dats(i).smooth.threshFreqC] = wgbss_smooth(dats(i).freqC(dats(i).rf),dats(i).pos(dats(i).rf),doplot);
    if(size(dats(i).smooth.freqC,1)< size(dats(i).smooth.freqC,2))
        dats(i).smooth.freqC = dats(i).smooth.freqC';
    end
    oldval = doSaveSmoothing;
    if doSaveSmoothing
        doSaveSmoothing = false;
    end
    [dats(i).smooth.nread(dats(i).rf),dats(i).smooth.threshNread] = wgbss_smooth(dats(i).nread(dats(i).rf),dats(i).pos(dats(i).rf),0);
    doSaveSmoothing = oldval;
    if(size(dats(i).smooth.nread,1)< size(dats(i).smooth.nread,2))
        dats(i).smooth.nread = dats(i).smooth.nread';
    end
end
