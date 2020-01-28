% Konstruktor der Klasse microarray-Experiment
% 
% E = meExperiment(name, file, chip, files)
%
%   name    Name des Projekts
% 
%   file = 'affy': filedata wird auf 
%       files.data  oder 'data_rma.csv'
%       files.calls oder 'data_calls.csv'
%   file = 'cdna': filedata wird auf 'Gene_Mean_Log_Ratio.xls' gesetzt.
%   file = []: filedata wird auf [pwd,'Gene_Mean_Log_Ratio.xls'] gesetzt.
%   file gar nicht �bergeben: filedata wird auf '@maExperiment/private/Gene_Mean_Log_Ratio.xls' gesetzt.
%   
%   chip  Bei affy: chip-Bezeichnung z.B. 'hgu133a' oder 'mouse4302'
% 
%   files.data
%   files.calls
% 
%   Beispiele:
%   E = maExperiment('Name',[]); 
%   E = maExperiment('Pagenstecher','pagenstecher__.rel')
%   E = maExperiment([ngene, na])
%   E = maExperiment(datenmatrix)
%   E = maExperiment(struct)
% 
%   E = maExperiment('name','*.gpr')
%   E = maExperiment('name','affy','mouse4302',struct('data','data_gcrma.csv','calls','data_calls.csv'))


function E = maExperiment(name, file, chip, files)    
global DataPath
if(~exist('files','var'))
    files = struct;
end

if(~exist('name','var'))
    pfad =which('maExperiment');
    load([pfad(1:end-14),filesep,'private',filesep,'Estruct.mat'])
    E = maExperiment(Estruct);
elseif(isstruct(name)) %KOnvertierung von Struct zu Klasse
    % Es wird schon eine Structur �bergeben, die dann in die Klasse
    % umgewandelt wird.
    a = name.arrays;
    for i=1:length(a)
        if(~isa(a{i},'microarray'))
            a{i} = microarray(a{i});
        end
    end
    name.arrays = a;

    Estruct = maStruct(name);
    E = class(Estruct,'maExperiment');   
%     E = class(name,'maExperiment');
elseif(isnumeric(name)) % Falls Datenmatrix �bergeben wurde oder Gr��e einer Matrix (f�r Zufallsdaten)
    if(sum(abs(size(name)-[1,2]))==0) %% Zufallszahlen
        for i=1:name(2)
            Estruct.arrays{i} = microarray(randn(name(1),1),['Chip',num2str(i)]);
        end
    else
        for i=1:size(name,2)
            Estruct.arrays{i} = microarray(name(:,i),['Chip',num2str(i)]);
        end
    end        
    Estruct = maStruct(Estruct);
    E = maExperiment(Estruct);
else
	DataPath = ['d:',filesep,'clemens',filesep,'microarray',filesep,'r',filesep];
	
	if(isempty(name)) % falls []
		name = '';
	end
    
	if(~exist('file','var'))
        file = [strrep(which('maExperiment'),'maExperiment.m','private\'),'Gene_Mean_Log_Ratio.xls'];
	elseif(isempty(file))
        file = [pwd,'\Gene_Mean_Log_Ratio.xls'];
	end
    disp(['file=',file])

    E = maStruct;
	E.annotation    = [];
    E.arrays        = cell(0);
    E.container.DreiStrichIDs = cell(0);
    E.container.FunfStrichIDs = cell(0);
	E.chromloc      = [];
    E.chromosome    = cell(0);
    E.clonename     = cell(0);
	E.genegroups    = cell(0);
    E.go            = cell(0);
    E.goIDs         = cell(0);
    E.goNames       = cell(0);
    E.groups        = cell(0);
    E.hybrname      = cell(0);
    E.IDs           = cell(0);
    E.llid          = cell(0);
    E.name          = '';
    E.namen         = cell(0);
    E.present       = [];
    E.symbol        = cell(0);
    E.unigeneID     = cell(0);    
    
	if(strcmp(file,'affy')==1)
        if(~isfield(files,'data'))
            filedata = 'data_rma.csv';
        else
            filedata = files.data;
        end
        if(~isfield(files,'calls'))
            filecalls = 'data_calls.csv';
        else
            filecalls = files.calls;
        end
        disp(['filedata  = ',filedata])
        disp(['filecalls = ',filecalls]);
        if(~exist('chip','var') | isempty(chip))
            chip = 'hgu133a';
        end
        disp(['maExperiment.m: Chip = ',chip])
        iobion = 0;
    elseif(strcmp(file,'*.gpr')==1)
        disp('*.gpr werden gelesen...')
    else 
        iobion = 1;
	end
	if(isempty(strfind(file,filesep)) & strcmp(file,'affy')~=1)
        file = [pwd,filesep,file];
	end
		
	
	[dummy,dummy,ext] = fileparts(file);
	if(strcmp(ext,'.rel')==1)
        relfile = 1;
        iobion = 0;
        absfile= 0;
    elseif(strcmp(ext,'.abs')==1)
        absfile = 1;
        relfile = 0;
        iobion  = 0;
    else 
        absfile = 0;
        relfile = 0;
	end
	
	
	E.name  = name;
	
	disp(['Lese ',file]);
	
	if(strcmp(file,'affy')==1)
        %%%%%%%%%%%%%%%%%%%%
        %%% Erste m�glichkeit: Affydaten:
        [dat, hybrname, IDs, present] = ReadAffymetrixData(filedata, filecalls);
        E.hybrname  = hybrname;
        E.IDs       = IDs;
	elseif(relfile==1)
        [dat, hybrname, IDs, descr, snr, relError] = ReadExpressionistData(file);
        E.hybrname  = hybrname;
        E.IDs       = IDs;
        E.container.snr         = snr;
        E.container.description = descr;
        E.container.relerror    = relError;
    elseif(absfile==1)
        [dat, hybrname, IDs, descr, p] = ReadExpressionistAbsData(file);
        E.hybrname  = hybrname;
        E.IDs       = IDs;
        E.container.absentp     = p;
        E.container.description = descr;  
        E.present   = p<0.05;
	else
        %%%%%%%%%%%%%%%%%%%%
        %%% Zweite M�glichkeit: cDna-DAten:
        try 
            [dat, hybrname, IDs] = ReadGeneMeanLogRatio(file);
            E.hybrname  = hybrname;
            E.IDs       = IDs;
        catch
            disp(['Fehlermeldung cDNA: ',lasterr]);
            error('maExperiment.m: Einlesen der Daten fehlgeschlagen.')
        end
        
    end

    if(~isempty(strfind('"',IDs{1})))
        error('maExperiment.m: " in den IDs.')
    end
	
    for i=1:length(E.hybrname)
        E.arrays{i}  = microarray(dat(:,i), E.hybrname{i});
    end
	    
    if(strcmp(file((end-3):end),'affy')~=1 & ~exist('affy','var'))
        disp('Annotation ...')
        pfadAnno = which('Annotation.mat','-all');
        load(pfadAnno{1});
	
        zuordn = struct;
        
        fnGenes = {'unigeneID','namen','symbol','llid','chromosome','clonename','go','goIDs','goNames'};
        for ifg = 1:length(fnGenes)
            E.(fnGenes{ifg}) = cell(length(E.IDs),1);
        end
        for i=1:length(E.IDs)
            try
                E.unigeneID{i}  = e.unigeneID{e.zuordn.(E.IDs{i})};
                E.namen{i}      = e.namen{e.zuordn.(E.IDs{i})};
                E.symbol{i}     = e.symbol{e.zuordn.(E.IDs{i})};
                E.llid{i}       = e.llid{e.zuordn.(E.IDs{i})};
                E.chromosome{i} = e.chromosome{e.zuordn.(E.IDs{i})};
                E.clonename{i}  = e.clonename{e.zuordn.(E.IDs{i})};
                E.go{i}         = e.go{e.zuordn.(E.IDs{i})};
                E.goIDs{i}      = e.goIds{e.zuordn.(E.IDs{i})};
                E.goNames{i}    = e.goNames{e.zuordn.(E.IDs{i})};
                E.container.DreiStrichIDs{i}    = e.DreiStrichIDs{e.zuordn.(E.IDs{i})};
                E.container.FunfStrichIDs{i}    = e.FunfStrichIDs{e.zuordn.(E.IDs{i})};
	%             try
	%                 eval(['e.zuordn.',E.IDs{i},'=',num2str(e.zuordn.(E.IDs{i})),';']);
	%             end
            catch
	%             disp(E.IDs{i})
	%             disp(lasterr)
                E.unigeneID{i}  = 'Data not found';
                E.namen{i}      = 'Data not found';
                E.symbol{i}     = 'Data not found';
                E.llid{i}       = 'Data not found';
                E.chromosome{i} = 'Data not found';
                E.clonename{i}  = 'Data not found';
                E.go{i}         = 'Data not found';
                E.goIDs{i}      = 'Data not found';
                E.goNames{i}    = '';
            end
        end
        
        indleer = find(cellfun('isempty',E.symbol));
        for iii=1:length(indleer)
            E.symbol{indleer(iii)} = '';
        end
        
        affy = 0;
        disp('Annotation ... fertig.')
    else
        affy = 1; 
        iobion = 0;
        maStruct;
    end
	
	
	E.groups = FindGroups(E.hybrname,length(E.hybrname));
	
	if(exist('present','var'))
        E.present = present;
    end
    	    
    E.container.date = date;
    E.container.time = datestr(now);
    E.container.pwd  = pwd;
	
% 	save 'maExperiment'
	
	E = class(E,'maExperiment');
	
	% fieldnames(E)
	
	if(affy==1 | absfile==1)
        disp('Annotation ...')
%         try

warning('             E = GetAffyAnnotation(E,chip);  in maExperiment disabled')
            disp('Annotation ... fertig.')
%         catch
%             disp('Annotation fehlgeschlagen.')
%             disp(lasterr);
%         end
	end    
	
	try
        miame = LoadMiameSampleAnnotation;
        E = set(E,'miame',miame);
	catch
        disp('MIAME Sample Annotation konnte nicht geladen werden. Grund:')
        disp(lasterr);
	end
    
    if(iobion==1)
            % Intensit�ten:
        try
            
            int = LoadLexFiles(E);
            fn = fieldnames(int);
            cont = get(E,'container');
            for f = 1:length(fn)
                cont.(fn{f}) = int.(fn{f});
            end
            E = set(E,'container',cont);
            [m,a] = maPlot(E);
            E = set(E,'m',m);
            E = set(E,'a',a);
            
            E = FlagSaturatedSpots(E);
            FindInterestingFlags(E);
        catch
            disp(lasterr)
            warning('maExperiment.m: Keine Intensit�ten geladen!')
        end
    end
end

