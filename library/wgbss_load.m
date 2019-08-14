% dat = wgbss_load(file)
% dat = wgbss_load(dat)
% 
%   loads data struct from file and converts the fields to proper class.
%   
% maxpos    maximal number of position, additional ones will be cut away
% 
% onlyFirstSample   if true, only the first column of the data is kept.
%               This reduces comp. effort e.g. for smooting


function dat = wgbss_load(file_or_struct, maxNpos, onlyFirstSample)

if ischar(file_or_struct)
    file = file_or_struct;
    dat = wgbss_bed2mat(file);
elseif isstruct(file_or_struct)
    dat = file_or_struct;
else
    file_or_struct
    error('argument does not fit');
end

if ~exist('onlyFirstSample','var') || isempty(onlyFirstSample)
    onlyFirstSample = false;
end


if onlyFirstSample
    disp('Use only one column/sample ...')
    nx = length(dat.pos);
    
    fn = fieldnames(dat);
    for f=1:length(fn)
        if size(dat.(fn{f}),1) == nx 
            dat.(fn{f}) = dat.(fn{f})(:,1);
        end
    end
end


if ~exist('maxNpos','var') || isempty(maxNpos)
    maxNpos = Inf;
else
    maxNpos = min(length(dat.nread),maxNpos);
end

if ~isinf(maxNpos)
    fprintf('Keep max. %i positions.\n',maxNpos);
    
    nx = length(dat.pos);
    dat.pos = dat.pos(1:maxNpos,:);
    
    fn = fieldnames(dat);
    for f=1:length(fn)
        if size(dat.(fn{f}),1) == nx
            dat.(fn{f}) = dat.(fn{f})(1:maxNpos,:);
        end
    end
end

dat.nread = single(dat.nread);
dat.pos = single(dat.pos);

if ~isfield(dat,'freqC') && isfield(dat,'freqT');
    if max(dat.freqT)>1
        dat.freqC = 100-dat.freqT;
    else
        dat.freqC = 1-dat.freqT;
    end
end

if min(dat.freqC(:))< -90  & max(dat.freqC(:))>0  % falsch skaliert von -99 bis +1
    if max(dat.freqT(:))>1 
        disp('Recalculate dat.freqC')
        dat.freqC = 100 - dat.freqT;
    end
elseif min(dat.freqC(:))<0
    warning('min(dat.freqC(:))<0');
end

if ~isfield(dat,'smooth') || ~isstruct(dat.smooth)
    disp('Smoothing ...')
    dat = wgbss_SmoothDat(dat);
end

if isfield(dat,'smooth')
%     dat.smooth.pos = single(dat.smooth.pos);
    dat.smooth.nread = single(dat.smooth.nread);  % for saving memory
end

