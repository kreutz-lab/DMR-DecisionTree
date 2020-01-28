% simF = getSimulationFeatures(force)
% 
%   This function returns simulation features i.e. parameters used for
%   simulating realistic data.
% 
%   force       if false (default), then the workspace simulationFeatures.mat is
%               loaded. force=true enforces recalculation and overwriting
%               simulationFeatures.mat
% 
%   simF        table
%               

function simF = getSimulationFeatures(force)

if (~force || nargin==0)  && exist('simulationFeatures.mat','file')
    load simulationFeatures simF

else
    % here a the following struct is loaded (calculated in
    % Auswertung5_SimulationParameters.m):
    %     
    % arg = 
    %   struct with fields:
    % 
    %     context: {1×21 cell}
    %         arg: {1×21 cell}        
    
    load Auswertung1.mat ass truth names_meth names_context
    load ../../WGBSSuite4_KreutzSimulator_Sep17/Auswertung5_SimulationParameters.mat arg
    
    dat.F1 = [ass.F1]';
    dat.method = names_meth';
    dat.context = names_context';
    dat.file = cell(size(dat.context));
    [val, names] = args2predictor(arg);
    arg.context = strrep(arg.context,'_',' ');
    for i=1:length(dat.context)
        dat.file{i} = truth{i}.file;
        ii = strmatch(lower(dat.context{i}),lower(arg.context),'exact');
        if length(ii)==1
            for j=1:length(names)
                dat.(str2fieldname(names{j}))(i,1) = val(ii,j);
            end
        else
            dat.context{i}
            error('dat.context{i} not found.');
        end
    end
    
    if length(unique(cellfun(@length,struct2cell(dat))))~=1
        dat
        warning('length of the fields of dat seems to differ. Check!')
    end
    
    dat.filename = dat.file;
    for i=1:length(dat.filename)
        [~,dat.filename{i}] = fileparts(dat.filename{i});
    end
    
    simF = struct2table(dat);    
    simF.Properties.RowNames = strcat(simF.filename,'_',simF.method);    

    save simulationFeatures simF
end
