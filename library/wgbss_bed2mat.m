% suc = wgbss_bed2mat(file)

function dat = wgbss_bed2mat(file)

[pfad,filename] = fileparts(file); % file name without file extension
matfile = [pfad,filesep,filename,'.mat'];

fprintf('wgbss_bed2mat(%s) ...\n',file);
if exist(matfile,'file')
    fprintf('\n%s exists and is read instead to save computation time. \n If not intended, please move or delete %s.\n\n',matfile,matfile);
    load(matfile); % variable dat is loaded
    if ~exist('dat','var')
        error('Workspace %s does not contain variable ''dat''',matfile)
    end
else
    try
        fprintf('Proceeding %s ...', file)
        if ~exist(file,'file')
            dat = [];
            error('File %s does not exist.',file);
            return
        end
        dat = biseqReadData_Aethioma(file,[],false);  % all lines, no header

        [~,dat.rf] = sort(dat.pos);
        dat = wgbss_SmoothDat(dat);
        
        save(matfile,'-v7.3','dat');
        fprintf('\n');
    catch ERR
        rethrow(ERR)
    end
end


