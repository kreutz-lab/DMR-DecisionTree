cd E:\clemens\Repositories\DMR-DecisionTree
%%
addpath('library');

% file = 'Data/P.abies/plus65_CHH.txt';
% file = 'Data/P.abies/plus65_CHG.txt'; 
% file = 'Data/P.abies/plus65_CG.txt';
% file = 'Data/P.abies/minus65_CHH.txt';
% file = 'Data/P.abies/minus65_CHG.txt';
% file = 'Data/P.abies/minus65_CG.txt';


% file = 'Data/A.thaliana/A.th-CG_Chr2.txt';
% file = 'Data/A.thaliana/A.th-CHG_Chr2.txt';
% file = 'Data/A.thaliana/A.th-CHH_Chr2.txt';

% file = 'Data/Ae.arabicum/Ae-T-CHG_Scaffold65.txt';
% file = 'Data/Ae.arabicum/Ae-T-CHH_Scaffold65.txt';
% file = 'Data/Ae.arabicum/Ae-c20-CG_Scaffold65.txt';
% file = 'Data/Ae.arabicum/Ae-C-CG_Scaffold65.txt';
% file = 'Data/Ae.arabicum/Ae-C-CHG_Scaffold65.txt';
% file = 'Data/Ae.arabicum/Ae-C-CHH_Scaffold65.txt';

% file = 'Data/P.patens/PP-G_CHG.Chr27.txt';
% file = 'Data/P.patens/PP-G_CHH.Chr27.txt';
% file = 'Data/P.patens/PP-G_CpG.Chr27.txt';
% file = 'Data/P.patens/PP-R_CHG.Chr27.txt';
% file = 'Data/P.patens/PP-R_CHH.Chr27.txt';
% file = 'Data/P.patens/PP-R_CpG.Chr27.txt';

dat = wgbss_load(file);                 % reading data
dat = wgbss_CalculateDataFeatures(dat); % calculation of data attributes

[X,categ,datF] = dat2X(dat,tree);       % formatting the data attributes as required for the decision tree
load Decision_tree.mat tree methodnames % loading the decision tree as matlab variable
 
pred = predict(tree,X)
fprintf('Decision tree suggests %s.\n',methodnames{round(pred)});

%% Test, ob der tree macht, was im Paper steht:

% view(tree,'Mode','graph')
i1 = strmatch('Meth_meanNread',datF.Properties.VariableNames,'exact');
i2 = strmatch('sdLog10MethLengthPos',datF.Properties.VariableNames,'exact');
i3 = strmatch('Meth_sd',datF.Properties.VariableNames,'exact');

X([i1,i2,i3])

xtest = X;
xtest(i1) = 3;
xtest(i2) = 0.5;
xtest(i3) = 40;
xtest([i1,i2,i3])
pred = predict(tree,xtest)




