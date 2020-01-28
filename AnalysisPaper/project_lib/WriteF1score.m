% WriteF1score(file,ass2,colnames)
% 
%   This function writes the plotted results as table (as requested by
%   Nilay)

function WriteF1score(file,ass2,colnames)

fid = fopen(file,'w');

for i=1:size(ass2,1)
    fprintf(fid,'%s\t%s\t',['true diff for ',colnames{i}],['F1 for ',colnames{i}]);
end
fprintf(fid,'\n');

for j=1:size(ass2,2) % over all difflev
    for i=1:size(ass2,1) % over all methods
        if ~isempty(ass2{i,j})
            out = sprintf('%f\t%f\t',ass2{i,j}.difflev,ass2{i,j}.F1);
            fprintf(fid,'%s',strrep(out,'.',','));
        else
            fprintf(fid,'\t\t');
        end
    end
    fprintf(fid,'\n');
end

fclose(fid);

