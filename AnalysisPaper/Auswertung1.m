clear all
addpath('project_lib')

%% location of the DMR results
datapath = 'DMR_Results';
addpath(datapath);

%% location of the simulation data *.mat files
simpath = 'BenchmarkData\BenchmarkData';  % has to be downloaded from figshare https://doi.org/10.6084/m9.figshare.11619045

d = dir(datapath);
files = setdiff({d.name},{'.','..'});

files_truth = cell(size(files));
names = cell(size(files));
names_meth = cell(size(files));
names_context = cell(size(files));
clc
for i=1:length(files)
    % replacements in order to have matching filename
    files_truth{i} = strrep(files{i},'_M1','');
    files_truth{i} = strrep(files_truth{i},'C20P','C20');
    files_truth{i} = strrep(files_truth{i},'C25P','C25');
    files_truth{i} = strrep(files_truth{i},'T20P','T20');
    files_truth{i} = strrep(files_truth{i},'T25P','T25');
    files_truth{i} = strrep(files_truth{i},'_200bp','');
    
    ind = strfind(files_truth{i},'_');
    [~,~,ext] = fileparts(files_truth{i});
    if strcmp(ext,'.bed')
        files_truth{i} = [files_truth{i}((ind(end-1)+1):(end-4)),'_SimulationData.mat'];
    elseif strcmp(ext,'.bedgraph')
        files_truth{i} = [files_truth{i}((ind(end-1)+1):(end-9)),'_SimulationData.mat'];
    else
        error('Extension %s not known.',ext);
    end
    
    % check whether for each DMR result the respective simulation data set
    % is found:
    if exist([simpath,filesep,files_truth{i}],'file')
        files_truth{i} = [simpath,filesep,files_truth{i}];
        fprintf('%s found.\n',files_truth{i})
    else
        error('%s not found.\n',files_truth{i})
    end
    
    % generate nice names:
    tmp = strsplit(files{i},'_');
    tmp2 = strsplit(tmp{end},'.');
    names_meth{i} = tmp{1};
    names_context{i} = [tmp{end-1},' ',tmp2{1}];
    names{i} = [names_meth{i},' & ',names_context{i}];
end


%% Load DMRs and truth
disp('Loading DMR calling results:')
res = cell(size(files));
for i=1:length(files)
    fprintf('importdata %s (%i out of %i) ...\n',files{i},i,length(files));
    res{i}.file = files{i};
    tmp = importdata(res{i}.file);
    res{i}.DMRs = tmp.data;
end

for i=1:length(res)
    res{i}.DMRs = res{i}.DMRs(:,1:2); % only the first two columns are relevant
end

%%
disp('Loading the simulation truth:')
truth = cell(size(files_truth));
for i=1:length(files_truth)
    fprintf('loading *.mat %s (%i out of %i) ...\n',files_truth{i},i,length(files_truth));
    truth{i}.file = files_truth{i};
    tmp = load(truth{i}.file);
    truth{i}.isDiff = tmp.sim.dat.isDiff;
    truth{i}.pos = tmp.sim.dat.pos;
end

save Auswertung1_tmp truth res -v7.3  % workspace with intermediate variables

%% Match DMR positions with CpG positions:
% Predicted DMRs do not always exactly coincide with pos in truth.
% res{i}.start_pos and res{i}.end_pos coincide.
for i=1:length(res)
    if ~isempty(res{i})
        fprintf('.')
        res{i}.start_pos = NaN(size(res{i}.DMRs,1),1);
        res{i}.end_pos = NaN(size(res{i}.DMRs,1),1);
        
        for j=1:size(res{i}.DMRs,1)
            smaller = find(truth{i}.pos < res{i}.DMRs(j,1));
            if isempty(smaller)
                smaller = 0;
            end
            if (smaller(end)+1)>length(truth{i}.pos)
                res{i} = [];
                break
            end
            res{i}.start_pos(j) = truth{i}.pos(smaller(end)+1);
            
            larger = find(truth{i}.pos > res{i}.DMRs(j,2));
            if isempty(larger)
                larger = length(truth{i}.pos)+1;
            end
            res{i}.end_pos(j) = truth{i}.pos(larger(1)-1);
        end
    end
end
fprintf('\n');

indOk = find(~cellfun(@isempty,res));
indNotOk =  find(cellfun(@isempty,res));
if ~isempty(indNotOk)
    fprintf('Empty: %s\n',sprintf(' %i ',indNotOk));
end
truth = truth(indOk);
res = res(indOk);
files = files(indOk);
names_contextAll = names_context;
names_context = names_context(indOk);
names = names(indOk);
names_methAll = names_meth;
names_meth = names_meth(indOk);


lev_meth = levels(names_meth);
index_meth = NaN(size(names_meth));
for i=1:length(names_meth)
    index_meth(i) = strmatch(names_meth{i},lev_meth,'exact'); % required for ensuring same colors
end

lev_context = levels(lower(names_context));
index_context = NaN(size(names_context));
for i=1:length(names_context)
    index_context(i) = strmatch(lower(names_context{i}),lower(lev_context),'exact'); % required for ensuring same colors
end

save res2 res truth indOk index_context index_meth -v7.3



%% create predicted DMPs from DMRs
for i=1:length(truth)
    fprintf('.')
    truth{i}.est.startMinusEnd = zeros(size(truth{i}.pos));
    [~,ia] = intersect(truth{i}.pos,res{i}.start_pos);
    truth{i}.est.startMinusEnd(ia) = 1;
    
    [~,ia] = intersect(truth{i}.pos,res{i}.end_pos);
    truth{i}.est.startMinusEnd(ia) = truth{i}.est.startMinusEnd(ia)-1;
    
    truth{i}.est.isDiff = cumsum(truth{i}.est.startMinusEnd);
end
fprintf('\n');
%% Design matrix
meth = struct;
meth.DMRcate = ~cellfun(@isempty,regexp(files,'DMRcate'));
meth.defiant = ~cellfun(@isempty,regexp(files,'defiant'));
meth.bsmooth = ~cellfun(@isempty,regexp(files,'bsmooth'));
meth.methylkit = ~cellfun(@isempty,regexp(files,'methylkit'));
meth.methylsig = ~cellfun(@isempty,regexp(files,'methylsig'));
meth.metilene = ~cellfun(@isempty,regexp(files,'metilene'));
meth.moabs = ~cellfun(@isempty,regexp(files,'moabs'));
meth.methylscore = ~cellfun(@isempty,regexp(files,'methylscore'));

context = struct;
context.G129_CHG = ~cellfun(@isempty,regexp(files,'G129_CHG'));
context.G129_CHH = ~cellfun(@isempty,regexp(files,'G129_CHH'));
context.G129_CpG = ~cellfun(@isempty,regexp(files,'G129_CpG'));
context.R143_CHG = ~cellfun(@isempty,regexp(files,'R143_CHG'));
context.R143_CHH = ~cellfun(@isempty,regexp(files,'R143_CHH'));
context.R143_CpG = ~cellfun(@isempty,regexp(files,'R143_CpG'));
context.Plus65_CpG = ~cellfun(@isempty,regexp(lower(files),'plus65_cpg'));
context.Plus65_CHG = ~cellfun(@isempty,regexp(lower(files),'plus65_chg'));
context.Plus65_CHH = ~cellfun(@isempty,regexp(lower(files),'plus65_chh'));
context.Minus65_CpG = ~cellfun(@isempty,regexp(lower(files),'minus65_cpg'));
context.Minus65_CHG = ~cellfun(@isempty,regexp(lower(files),'minus65_chg'));
context.Minus65_CHH = ~cellfun(@isempty,regexp(lower(files),'minus65_chh'));
context.C20_CpG = ~cellfun(@isempty,regexp(lower(files),'c20_cpg'));
context.C20_CHG = ~cellfun(@isempty,regexp(lower(files),'c20_chg'));
context.C20_CHH = ~cellfun(@isempty,regexp(lower(files),'c20_chh'));
context.C25_CpG = ~cellfun(@isempty,regexp(lower(files),'c25_cpg'));
context.T20_CHG = ~cellfun(@isempty,regexp(lower(files),'t20_chg'));
context.T20_CHH = ~cellfun(@isempty,regexp(lower(files),'t20_chh'));
context.arabidopsis_CpG = ~cellfun(@isempty,regexp(lower(files),'arabidopsis_cpg'));
context.arabidopsis_CHG = ~cellfun(@isempty,regexp(lower(files),'arabidopsis_chg'));
context.arabidopsis_CHH = ~cellfun(@isempty,regexp(lower(files),'arabidopsis_chh'));


fn = fieldnames(meth);
Xmeth = zeros(length(meth.(fn{1})),length(fn));
for i=1:length(fn)
    Xmeth(:,i) = meth.(fn{i});
end

fn = fieldnames(context);
Xcontext = zeros(length(context.(fn{1})),length(fn));
for i=1:length(fn)
    Xcontext(:,i) = context.(fn{i})';
end

%% Reorder contexts
context2_names = {'CpG','CHH','CHG'};
fnc = fieldnames(context);
Xcontext2 = [];
for i=1:length(context2_names)
    ind = find(~cellfun(@isempty,strfind(fnc,context2_names{i})));    
    Xcontext2(:,i) = sum(Xcontext(:,ind),2);
end

%% Simple features:
clear CPs ass
for i=1:length(truth)
    x = truth{i}.isDiff~=0;
    y = truth{i}.est.isDiff;
    
    ind_analyzed = 1:length(truth{i}.pos);    
    % the range 50000 - 100000 has to be excluded (see email 28.9.18)
    if meth.bsmooth(i)==1 && ~isempty(strmatch('arabidopsis',names_context{i}))
        ind_analyzed = find(truth{i}.pos<50000 | truth{i}.pos>100000);
        names_context{i}
    end
    truth{i}.ind_analyzed = ind_analyzed;
    truth{i}.bool_analyzed = false(size(truth{i}.pos));
    truth{i}.bool_analyzed(truth{i}.ind_analyzed) = true;
    
    if i==1
        ass = assess_classification(x(ind_analyzed),y(ind_analyzed));
    else
        ass(i) = assess_classification(x(ind_analyzed),y(ind_analyzed));
    end
end

%% linear model:
refmeth = 'metilene'
fnmeth = fieldnames(meth);
iref = strmatch(refmeth,fnmeth,'exact');
iNotRef = setdiff(1:length(fnmeth),iref);

refcontext = 'R143_CpG';
fncon = fieldnames(context);
iref2 = strmatch(refcontext,fncon,'exact');
iNotRef2 = setdiff(1:length(fncon),iref2);

if isempty(iref) || isempty(iref2)
    error('isempty(iref) || isempty(iref2): refmeth or refcontext not found.')
end

xnames = [{'Reference'},fnmeth(iNotRef)',fncon(iNotRef2)'];
refname = [refmeth,' & ',refcontext];

[b,bint,r,rint,stats] = regress([ass.F1]',[ones(size(Xmeth,1),1),Xmeth(:,iNotRef),Xcontext(:,iNotRef2)]);
mdl = fitlm([ones(size(Xmeth,1),1),Xmeth(:,iNotRef),Xcontext(:,iNotRef2)],[ass.F1]','PredictorVars',xnames,'Intercept',false) % for p-values
s = struct(mdl);
tmp = table2array(s.Coefficients);
pvals = tmp(:,end);
siglabel = p2siglabel(pvals);
for i=1:length(siglabel)
    siglabel{i} = sprintf('%5s',deblank(siglabel{i}));
end
WriteTableBioinformatics('Multivariate.tex',table2array(s.Coefficients),{'estimate','SE','t-statistic','p-value'},ReplaceSampleNames(strrepPaper(s.CoefficientNames)));

[b2,bint2,r2,rint2,stats2] = regress([ass.F1]',[ones(size(Xmeth,1),1),Xmeth(:,iNotRef),Xcontext2]);


%%
close all
patch([0.5,1.5,1.5,.5],[-.6,-.6,1,1],ones(1,3)*0.8);
hold on
doplot = 1:length(b); 
doplot(strmatch('methylscore',xnames,'exact')) = [];
for i=1:length(doplot)
    if i==1
        plot([i,i],bint(doplot(i),:),'LineWidth',2,'Color',ones(1,3)*0)
        plot([i],b(doplot(i)),'o','LineWidth',2,'Color',ones(1,3)*0)        
    else
        plot([i,i],bint(doplot(i),:),'k','LineWidth',2)
        plot([i],b(doplot(i)),'ko','LineWidth',2)
    end
end
xlim([.5,length(b)+.1])
abplot(0,0)
xticklabel_rotate(1:length(doplot),90,strcat(ReplaceSampleNames(str2label(strrepPaper(xnames(doplot)))),siglabel(doplot)'));
ylabel('Impact on F1-score, 95% CIs');
title(str2label(['Reference analysis: ',ReplaceSampleNames(strrep(strrepPaper(refname),'-&-',' & '))]))
ylim([-.6,1])
print -dpng MultivariateStatisticsB

%% depending on the magnitude of regulation:
clear ass2
for i=1:length(truth)
    difflev = unique(truth{i}.isDiff(truth{i}.bool_analyzed));
    fprintf('.');
    
    if length(difflev)>100
        difflev = linspace(min(difflev),max(difflev),100); % reduce comp effort
    end
    
    for j=1:length(difflev)
        ind = find(truth{i}.isDiff(truth{i}.bool_analyzed)<=(difflev(j)+1e-5));
        
        x = truth{i}.isDiff(truth{i}.bool_analyzed);
        x = x(ind)~=0;
        
        y = truth{i}.est.isDiff(truth{i}.bool_analyzed);
        y = y(ind);
        
        tmp = assess_classification(x,y,true); % fast option, only F1
        tmp.difflev = difflev(j);
        ass2{i,j} = tmp;
    end
    
end
fprintf('\n');

WriteF1score('F1-scores.xls',ass2,strcat(names_context,'_',names_meth))

% save tmp -v7.3
%%
ass2x = NaN(size(ass2));
ass2y = NaN(size(ass2));
for i=1:size(ass2,1)
    for j=1:size(ass2,2)
        if ~isempty(ass2{i,j})
            ass2x(i,j) = ass2{i,j}.difflev;
        end
        if ~isempty(ass2{i,j})
            ass2y(i,j) = ass2{i,j}.F1;
        end
    end
end

%%
close all
fncon = fieldnames(context);
nx = ceil(sqrt(length(fncon)));
ny = ceil(length(fncon)/nx);
load col % load colors
close all
phase_diff = 0.25;  % x-axis has to be scaled with phase_diff
for i=1:length(fncon)
    figure(1)
    subplot(ny,nx,i)

    iy = find(context.(fncon{i})==1);
    if ~isempty(iy)
        iys(i) = iy(1);
        doplot = 1:length(iy);
        doplot(strmatch('methylscore',names_meth(iy),'exact')) = []; 
        for c=doplot 
            plot(phase_diff*ass2x(iy(c),:)',ass2y(iy(c),:)','.-','Color',col(index_meth(iy(c)),:));
            hold on
        end
        if i==1
            hleg = legend(names_meth(iy(doplot)),'Location','NorthWest');
        else
%             legend(names_meth(iy),'Location','NorthWest');     % can be uncommented for checking the colors       
        end
        title(strrepPaper(names_context{iys(i)}));
        
        figure(2)
        for c=doplot%1:length(iy)
            plot(phase_diff*ass2x(iy(c),:)',ass2y(iy(c),:)','.-','Color',col(index_meth(iy(c)),:),'LineWidth',2);
            hold on
        end
        legend(names_meth(iy(doplot)),'Location','NorthWest');     % can be uncommented for checking the colors
        title(strrepPaper(names_context{iys(i)}));
        ylabel('F1-score');
        xlabel('true diff. methylation');
        set(gca,'FontSize',16,'LineWidth',2);
        PrintToPng(gcf,names_context{iys(i)});
        close(2)
                
    else
        set(gca,'XTick',[],'YTick',[]);
        title([str2label(fncon{i}),': not yet analyzed'],'FontSize',8);
    end
end
suplabel('F1-score','y');
suplabel('true differential methylation','x');
set(gcf,'Position',1e3*[0.0081    0.1050    1.004    0.6529])
set(hleg,'Position',[0.6278    0.1233    0.1058    0.1480])
saveas(gcf,'DependencyOnDiffMeth')
print -dpng DependencyOnDiffMeth

%%
xinter = linspace(0,max(ass2x(:)),size(ass2y,2));
ass2yInter = NaN(size(ass2y));
for i=1:size(ass2x,1)
    notnan = find(~isnan(ass2x(i,:)) & ~isnan(ass2y(i,:)));
    ass2yInter(i,:) = interp1(ass2x(i,notnan),ass2y(i,notnan),xinter,'linear','extrap');
end

figure
hold on
for m=doplot 
    ind = find(index_meth==m);
    plot(xinter*phase_diff,mean(ass2yInter(ind,:),1)','.-','Color',col(m,:),'LineWidth',2);
end
legend(ReplaceSampleNames(names_meth(iy(doplot))),'Location','NorthWest');     % can be uncommented for checking the colors
ylabel('F1-score');
ylim([0,0.8])
xlabel('true differential methylation level');
set(gca,'FontSize',16,'LineWidth',2);
PrintToPng(gcf,'DependencyOnDiffMethAll');

%%
for plotm=1:length(lev_meth)
    figure
    hold on
    h = 0;
    for m=plotm
        ind = find(index_meth==m);
        for i=1:length(ind)
            h(m) = plot(xinter*phase_diff,ass2yInter(ind(i),:)','.-','Color',col(m,:),'LineWidth',2);
        end
    end
    legend(h(plotm),lev_meth(plotm),'Location','NorthWest');     % can be uncommented for checking the colors
    ylabel('F1-score');
    ylim([0,1])
    xlabel('true differential methylation level');
    set(gca,'FontSize',16,'LineWidth',2);
    PrintToPng(gcf,['DependencyOnDiffMeth_',lev_meth{m}]);
end

%% Rank image
close all
[rang,F1,rangSorted] = ImageBest(ReplaceSampleNames(names_meth),strrepPaper(names_context),ass2x,ass2y,[1,2]);

%%
save Auswertung1 -v7.3

