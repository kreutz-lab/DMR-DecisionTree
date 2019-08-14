# DMR-DecisionTree
Data-based prediction of the best performing DMR algorithm

## Code example
The decision tree requires data attributes like the average methylation levels in presumably methylated/unmethylated regions. 
This features can be calulated by the code provided in the repository as shown in the following example:

```
addpath('library');          % adding the folder to Matlab search path

dat = wgbss_load(file);      % reading data
dat = wgbss_CalculateDataFeatures(dat); % calculation of data attributes

X = dat2X(dat,tree);         % formatting the data attributes as required for the decision tree
load Decision_tree.mat tree methodnames % loading the decision tree and method names as matlab variable

pred = predict(tree,X);      % performing predictions
fprintf('Decision tree suggests %s.\n',methodnames{round(pred)});  % printing the outcome
```

## Data format
Our code requires the following data format:
````
Scaffold_65	16310	+	2	100.000
Scaffold_65	16332	+	1	100.000
Scaffold_65	16357	+	2	100.000
Scaffold_65	16379	+	2	100.000
Scaffold_65	16466	+	3	0.000
Scaffold_65	16471	+	2	0.000
Scaffold_65	16908	-	1	0.000
Scaffold_65	17545	+	1	0.000
Scaffold_65	17588	+	3	0.000
Scaffold_65	17626	+	2	50.000
Scaffold_65	17647	+	3	33.333
Scaffold_65	17673	+	5	20.000
Scaffold_65	17675	-	2	0.000
Scaffold_65	17716	+	6	0.000
Scaffold_65	17717	-	1	0.000
```
Description of the columns:
1) Chrosome (or scaffold) name
2) Position [bp]
3) Strand
4) Reads [counts]
5) Methylation level [%]

