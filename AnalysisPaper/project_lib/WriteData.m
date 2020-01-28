% Schreibt das maExperiment E in eine Datei mit Namen file.
% 
% status = WriteData(E,file)

function status = WriteData(E,file,decsep)
if(~exist('decsep','var') | isempty(decsep))
    decsep = '.';
end


IDs        = get(E,'IDs');
unigeneID  = get(E,'unigeneID');
if(isempty(unigeneID))
    unigeneID = empty_cellstring(size(IDs,1),1,'');
end
namen      = get(E,'namen');
if(isempty(namen))
    namen = empty_cellstring(size(IDs,1),1);
end
symbol     = get(E,'symbol');
if(isempty(symbol))
    symbol = empty_cellstring(size(IDs,1),1);
end
llid       = get(E,'llid');
chromosome = get(E,'chromosome');
go         = get(E,'go');

data       = get(E,'data');
hybrname   = get(E,'hybrname');

% titel = 'GeneBank Accession ID\tUnigene ClusterID\tGenename\tsymbol\tLLID\tchromosome\tGO';
titel = 'GeneBank Accession ID\tUnigene ClusterID\tGenename\tsymbol';
for i=1:length(hybrname)
    titel = [titel,'\t',hybrname{i}];
end
titel = [titel,'\n'];

fid = fopen(file,'w');
fprintf(fid,titel);
for ig=1:length(IDs)
    fprintf(fid,'%s',cell2mat(IDs(ig)));
    fprintf(fid,'\t%s',cell2mat(unigeneID(ig)));
    fprintf(fid,'\t%s',cell2mat(namen(ig)));
    fprintf(fid,'\t%s',cell2mat(symbol(ig)));
%     fprintf(fid,'\t%s',cell2mat(llid(ig)));
%     fprintf(fid,'\t%s',cell2mat(chromosome(ig)));
%     fprintf(fid,'\t%s',cell2mat(go(ig)));
    for ih = 1:length(hybrname)
        fprintf(fid,'\t%s',strrep(sprintf('%f',data(ig,ih)),'.',decsep));
    end
    fprintf(fid,'\n');
end

status = fclose(fid);




