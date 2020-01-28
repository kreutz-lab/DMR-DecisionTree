% status = WriteTableBioinformatics(file, data, colnames)
% status = WriteTableBioinformatics(file, data, colnames, rownames, upperleftname, nrboldbest, align, label, landscape)
% 
%   Schreibt den Latex-Code für eine Tabelle im Bioinformatics-Style
% 
% 
% status = WriteTableBioinformatics(file, data, colnames, rownames, ...
%           upperleftname, nrboldbest, align, label, landscape)
% 
%   struct  Die Felder müssen alles Vektoren der gleichen Länge sein.
%           Die Feldnamen werden als Spaltenbezeichnungen verwendet.
%           Existiert ein Feld "rownames", so wird es als Bezeichnung der
%           Zeilen verwendet.
%   rownames    Im Falle eines Structs werden die Unterstriche '_' werden
%               ersetzt durch Leerzeichen. 
%   nrboldbest  Anzahl der besten Samples, die durch Fettdruck
%               hervorgehoben werden.
%               Negative Zahlen: Die kleinsten.
%               Positive Zahlen: Die groessten.
%               Array der länge size(data,2)
%   align       Die Ausrichtung der Spalten, 
%               String der länge size(data,2)+1
%   label       Die Tabellen-Unterschrift evtl. mit dem Label
%   landscape   Default: 0
%               Bei landscape=1 wird die Tabelle ins Querformat gedreht.

function status = WriteTableBioinformatics(file, data, colnames, rownames, upperleftname, nrboldbest, align, label, landscape)
if(~exist('upperleftname','var') | isempty(upperleftname))
    upperleftname = ' ';
end
if(~exist('nrboldbest','var') | isempty(nrboldbest))
    nrboldbest = [];
end
if(~exist('align','var') | isempty(align))
    align = 'l';
    for i=1:size(data,2)
        align = [align,'l'];
    end
end
if(~exist('label','var') | isempty(label))
    label = '';
end
if(~exist('landscape','var') | isempty(landscape))
    landscape = 0;
end

%% Falls Übergabe als Struct:
if(isstruct(data))
    s = data;
    clear data;
    if(isfield(s,'rownames'))
        rownames = s.rownames;
        s = rmfield(s,'rownames');
    end
    colnames = fieldnames(s);
    data = NaN*ones(length(s.(colnames{1})),length(colnames));
    for i=1:length(colnames)
        data(:,i) = s.(colnames{i})(1:end)';
    end    
end

if(~iscell(data))
    data_tmp = cell(size(data));
    for i1=1:size(data,1)
        for i2=1:size(data,2)
            data_tmp{i1,i2} = data(i1,i2);
        end
    end
    data = data_tmp;
end


%% Schreiben:
fid = fopen(file,'w');
% fprintf(fid,'%s\n','\documentclass{bioinfo}');
% fprintf(fid,'%s\n','\copyrightyear{2006}');
% fprintf(fid,'%s\n','\pubyear{2006}');
% fprintf(fid,'%s\n','\usepackage{lscape}');
% fprintf(fid,'%s\n','\begin{document}');
if(landscape==1)
    fprintf(fid,'%s\n','\begin{landscape}');
    fprintf(fid,'%s\n','\begin{table}[!t]');   
else
    fprintf(fid,'%s\n','\begin{table*}[!t]');   
end

fprintf(fid,'%s%s%s\n','\processtable{',label,'}{');
fprintf(fid,'\t%s \n',['\begin{tabular}{',align,'}\toprule']);   

fprintf(fid,'%s\t & ',upperleftname);
for i=1:length(colnames)
    fprintf(fid,'%s %s\t',' ',colnames{i});
    if(i<length(colnames))
        fprintf(fid,' & ');
    else
        fprintf(fid,' %s\n ','\\\midrule');
    end
end

%% die besten ermitteln:
dobold = zeros(size(data));
for i=1:size(data,2)
    tmp = [data{:,i}];
    if(isnumeric(tmp) & ~isempty(nrboldbest))
        if(nrboldbest(i)<0)
            [tmpsort,rf] = sort(tmp);
            dobold(rf(1:abs(nrboldbest(i))),i) = 1;
        elseif(nrboldbest(i)>0)
            [tmpsort,rf] = sort(-tmp);
            dobold(rf(1:abs(nrboldbest(i))),i) = 1;            
        end
    end
end

%% Die Daten schreiben
for ig=1:size(data,1)
    fprintf(fid,'%s\t & ',rownames{ig});
    for ih = 1:size(data,2)
        if(dobold(ig,ih)==1)
            bfstr1 = '\bf{';
            bfstr2 = '}';
        else
            bfstr1 = '';
            bfstr2 = '';
        end
        if(isa(data,'cell') & isnumeric(data{ig,ih}))            
            fprintf(fid,'\t%s%s%s',bfstr1,strrep(dprintf(data{ig,ih}),'NaN','-'),bfstr2);
        elseif(isa(data,'cell'))
            fprintf(fid,'\t%s%s%s',bfstr1,strrep(data{ig,ih},'NaN','-'),bfstr2);
        else
            fprintf(fid,'\t%s%s%s',bfstr1,strrep(dprintf(data(ig,ih)),'NaN','-'),bfstr2);
        end
        if(ih<size(data,2))
            fprintf(fid,' & ');
        end
    end
    if(ig<size(data,1))
        fprintf(fid,'%s\n','\\');
    else
        fprintf(fid,'%s\n','\\\botrule');
    end
end

fprintf(fid,'%s\n','\end{tabular}}{}');   
if(landscape==1)
    fprintf(fid,'%s\n','\end{table}');   
    fprintf(fid,'%s\n','\end{landscape}');  
else
    fprintf(fid,'%s\n','\end{table*}');   
end
% fprintf(fid,'%s\n','\end{document}');  

status = fclose(fid);
