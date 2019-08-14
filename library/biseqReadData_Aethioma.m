function data = biseqReadData_Aethioma(file,nlines,hasheader,isParabies)
if ~exist('hasheader','var')
    hasheader = true;
end
if ~exist('isParabies','var') || isempty(isParabies)
    isParabies = false;
end
if ~exist('nlines','var') || isempty(nlines)
    nlines = wc(file);
end

[~,~,EXT] = fileparts(file);
EXT = strtrim(EXT);

if strcmp(EXT,'.txt')==1
    data = biseqReadData_txt_Aethioma(file,nlines,hasheader,isParabies);
elseif strcmp(EXT,'.bed')==1
    data = biseqReadData_bed(file,nlines,hasheader);
else
    EXT
    error('biseqReadData.m: reader only implemented for *.txt and *.bed');
end


function data = biseqReadData_bed(file,nlines,hasheader)
data.chr = cell(nlines,1);
data.pos = NaN(nlines,1);
data.nread = NaN(nlines,1);
data.strand = cell(nlines,1);
data.freqT = NaN(nlines,1);

fid = fopen(file,'r');
if(hasheader)
    textscan(fid, '%q %q %q %q %q %q %q\n',1);
end

i=0;
C = NaN;
while ~isempty(C)
    i=i+1;
    C = textscan(fid, '%q %f %f %f\n',1);
    
    data.chr{i} = C{1}{1};
%     data.pos(i) = C{2};
    data.pos(i) = C{3};
    data.strand{i} = 'F';  % no info provided, use default "F"
    data.nread(i) = NaN;
    data.freqT(i) = C{4};
    if(i>=nlines)
        break;
    end
end
fclose(fid);

data = biseqDataProperties(data);
data.file = file;
[~,data.sample]=fileparts(data.file);

data.CSpos = strcat(data.chr,'_',data.strand,'_',num2str(data.pos));




function data = biseqReadData_txt(file,nlines,hasheader)
data.chr = cell(nlines,1);
data.pos = NaN(nlines,1);
data.nread = NaN(nlines,1);
data.strand = cell(nlines,1);
data.freqT = NaN(nlines,1);
% disp(file)

fid = fopen(file,'r');
if(hasheader)
    textscan(fid, '%q %q %q %q %q %q %q\n',1);
end

i=0;
C = NaN;
while ~isempty(C)
    i=i+1;
    C = textscan(fid, '%q %q %f %q %f %f %f\n',1);

    cl = cellfun(@length,C);
    if(sum(cl([2:5,7])==0)==0) % wenn was drin steht
        try
            data.chr{i} = C{2}{1};
            data.pos(i) = C{3};
            data.strand{i} = C{4}{1};
            data.nread(i) = C{5};
            data.freqT(i) = C{7};
        catch
            fprintf('Problem with file %s: Line %i \n',file,i);
            C
            celldisp(C)
        end
    end
    
    if(i>=nlines)
        break;
    end
end
fclose(fid);

data = biseqDataProperties(data);
data.file = file;
[~,data.sample]=fileparts(data.file);

data.CSpos = strcat(data.chr,'_',data.strand,'_',num2str(data.pos));



function data = biseqReadData_txt_Aethioma(file,nlines,hasheader,isParabies)
data.chr = cell(nlines,1);
data.pos = NaN(nlines,1);
data.nread = NaN(nlines,1);
data.strand = cell(nlines,1);
data.freqT = NaN(nlines,1);
disp(file)

fid = fopen(file,'r');
if(hasheader)
    textscan(fid, '%q %q %q %q %q\n',1);
end

i=0;
C = NaN;
raus = ones(1,nlines);
while ~isempty(C)
    i=i+1;
    C = textscan(fid, '%q %f %q %f %f\n',1);

    cl = cellfun(@length,C);
    if(sum(cl([1:5])==0)==0) % wenn was drin steht
        try
            data.chr{i} = C{1}{1};
            data.pos(i) = C{2};
            data.strand{i} = C{3}{1};
            data.nread(i) = C{4}+C{5};  
            data.freqT(i) = 100*C{5}/(C{5}+C{4});
            raus(i) = 0;
        catch ERR
            fprintf('Problem with file %s: Line %i \n',file,i);
            C
            celldisp(C)
            rethrow(ERR)
        end
    end
    
    if(i>=nlines)
        break;
    end
end
fclose(fid);

fn = fieldnames(data);
for f=1:length(fn)
    data.(fn{f}) = data.(fn{f})(~raus);
end

if isParabies  % beinhaltet mehrer Chr, positionen müssen geschiftet werden
    dpos = [0;diff(data.pos)];
    q95 = quantile(dpos(dpos>0),.95); % zwischen Chr soll immer Abstand 95% quantile sein
    dpos(dpos>0) = 0;
    dpos(dpos<0) = -dpos(dpos<0)+q95;
    cdpos = cumsum(dpos);
    data.pos = data.pos+cdpos;
end


data = biseqDataProperties(data);
data = biseqMergeStrands(data);

data.file = file;
[~,data.sample]=fileparts(data.file);

data.CSpos = strcat(data.chr,'_',data.strand,'_',num2str(data.pos));

