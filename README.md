# DMR-DecisionTree
Data-based prediction of the best performing DMR algorithm

## Code example
The decision tree requires data attributes like the average methylation levels in presumably methylated/unmethylated regions. 
This features can be calulated by the code provided in the repository as shown in the following example:

```
addpath('library');          % adding the folder to Matlab search path

dat = wgbss_load(file);      % reading data
dat = wgbss_CalculateDataFeatures(dat); % calculation of data attributes

X = dat2X(dat);                         % formatting the data attributes as required for the decision tree
load Decision_tree.mat tree methodnames % loading the decision tree and method names as matlab variable

pred = predict(tree,X);      % performing predictions
fprintf('Decision tree suggests %s.\n',methodnames{round(pred)});  % printing the outcome
```
