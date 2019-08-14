% [X,categ] = dat2X(dat,tree)
% 
% This function calculates predictors from the data set
% 
%   X   predictors as table (i.e. design matrix format as in the original analysis)
% 
% Example:
% dat = wgbss_load('Data/A.thaliana/A.th-CG_Chr2.txt');
% dat = wgbss_CalculateDataFeatures(dat);
% X = dat2X(dat);

function [X,categ,datF] = dat2X(dat,tree)
s = struct;
s.val = wgbss_optimize_dat2res(dat);
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

    datF = struct;
    for i=1:length(s.names)
        datF.(s.names{i}) = s.val(:,i);
    end
    datF = struct2table(datF);
%     datF = [datF,table({file})];
%     datF.Properties.VariableNames{end} = 'file';
    datF = [datF,table({'NA'})];
    datF.Properties.VariableNames{end} = 'method';
    
    fn = setdiff(datF.Properties.VariableNames,{'F1','context','file','filename','target1','target2'});
X = [];
for i=1:length(fn)
    if isnumeric(datF.(fn{i}))
        X(:,i) = datF.(fn{i});
    else
        [~,levnr] = levels(datF.(fn{i}));
        X(:,i) = levnr;
    end
end

%% ensure that the predictors calculated in this function have the same order as in the decision tree:
warning off
s = struct(tree);
warning on

[~,ia,ib]=intersect(datF.Properties.VariableNames,s.PredictorNames);
if sum(diff(ib)<1)>0
    error('PredictorNames in the decision tree are not sorted. This is not a gernal issue but is assumed at this point.')
end
if length(ib)<length(datF.Properties.VariableNames)
    setdiff(s.PredictorNames,datF.Properties.VariableNames)
    error('Some predictors are missing. This is only a problem, if they are used by the tree. Fill the missing ones.');
end
datF = datF(:,ia);
X = X(:,ia);


%%
fn = setdiff(datF.Properties.VariableNames,{'F1','context','file','filename','target1','target2'});
categ = []; % whcih vars are categorical
for i=1:length(fn)
    if ~isnumeric(datF.(fn{i}))
        categ = [categ,i];
    end
end
