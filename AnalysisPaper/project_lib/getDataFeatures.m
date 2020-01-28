% datF = getDataFeatures(force)
% 
%   This function returns data features as used as residuals when making
%   the simulation data as "similar" as possible to real data
% 
% 
%   force       if false (default), then the workspace dataFeatures.mat is
%               loaded. force=true enforces recalculation and overwriting
%               dataFeatures.mat
% 
%   datF       table

function datF = getDataFeatures(force)

if (~force || nargin==0)  && exist('dataFeatures.mat','file')
    load dataFeatures 

else

    addpath('../../WGBSSuite4_KreutzSimulator_Sep17/project_lib/');
    load Auswertung1 files_truth names_methAll
    
    s = struct;
    s.val = [];
    for f=1:length(files_truth)
        fprintf('%i ',f);
        tmp = load(files_truth{f});
        tmp.sim.dat=wgbss_SmoothDat(tmp.sim.dat);
        tmp.sim.dat=wgbss_CalculateDataFeatures(tmp.sim.dat);
        
        s.val(f,:) = wgbss_optimize_dat2res(tmp.sim.dat);
        
        if rem(f,10)==0
            fprintf(' out of %i done.\n',length(files_truth));
        end
    end
    fprintf('\n');
    
    s.names = {'Meth_mean', 'NonMeth_mean', 'Meth_sd', 'NonMeth_sd', ...
        'Meth_meanNread', 'Meth_sdNread', ...
        'NonMeth_meanNread', 'NonMeth_sdNread', ...
        'meanLog10DiffPos', 'sdLog10DiffPos', ...
        'meanLog10MethLengthPos', 'sdLog10MethLengthPos', ...
        'meanLog10NonMethLengthPos', 'sdLog10NonMethLengthPos', ...
        'smooth_threshFreqC1',...
        'smooth_threshFreqC2',...
        'smooth_threshNread1',...
        'smooth_threshNread2'};
    
    for i=1:length(s.names)
        datF.(s.names{i}) = s.val(:,i);
    end
    
    datF = struct2table(datF);
    datF = [datF,table(files_truth')];
    datF.Properties.VariableNames{end} = 'file';
    
    datF.filename = datF.file;
    for i=1:length(datF.filename)
        [~,datF.filename{i}] = fileparts(datF.filename{i});
    end
    datF.Properties.RowNames = strcat(datF.filename,'_',names_methAll');
    
    save dataFeatures datF 
end


