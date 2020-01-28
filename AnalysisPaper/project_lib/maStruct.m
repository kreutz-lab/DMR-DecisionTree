% E = maStruct(Ein)
% 
%    Erzeugt einen Struct, der die Felder für die Klasse maExperiment hat.
% 
%   Ein  Kann struct sein, der schon bestimmte Felder besitzt.
%        Reihenfolge muss noch nicht stimmen.
% 
% Bsp:
    % e.arrays{1} = microarray(randn(100,1));
    % Estruct = maStruct(e);
    % E = maExperiment(Estruct);

function E = maStruct(Ein)
if(~exist('Ein','var'))
    Ein = struct;
end

if(~isfield(Ein,'name'))
    E.name = '';
else
    E.name = Ein.name;
end
if(~isfield(Ein,'container'))
    E.container = struct;
else
    E.container = Ein.container;
end

if(~isfield(Ein,'hybrname'))
    E.hybrname = [];
else
    E.hybrname = Ein.hybrname;
end
if(~isfield(Ein,'IDs'))
    E.IDs = [];
else
    E.IDs = Ein.IDs;
end

if(~isfield(Ein,'arrays'))
    E.arrays{1} = microarray([]);
else
    E.arrays = Ein.arrays;
end


if(~isfield(Ein,'unigeneID'))
    E.unigeneID = [];
else
    E.unigeneID = Ein.unigeneID;
end
if(~isfield(Ein,'namen'))
    E.namen = [];
else
    E.namen = Ein.namen;
end
if(~isfield(Ein,'symbol'))
    E.symbol = [];
else
    E.symbol = Ein.symbol;
end
if(~isfield(Ein,'llid'))
    E.llid = [];
else
    E.llid = Ein.llid;
end
if(~isfield(Ein,'chromosome'))
    E.chromosome = [];
else
    E.chromosome = Ein.chromosome;
end
if(~isfield(Ein,'clonename'))
    E.clonename = [];
else
    E.clonename = Ein.clonename;
end
if(~isfield(Ein,'go'))
    E.go = [];
else
    E.go = Ein.go;
end
if(~isfield(Ein,'goIDs'))
    E.goIDs = [];
else
    E.goIDs = Ein.goIDs;
end
if(~isfield(Ein,'goNames'))
    E.goNames = [];
else
    E.goNames = Ein.goNames;
end

if(~isfield(Ein,'groups'))
    E.groups = cell(0);
else
    E.groups = Ein.groups;
end
if(~isfield(Ein,'present'))
    E.present = [];
else
    E.present = Ein.present;
end
if(~isfield(Ein,'chromloc'))
    E.chromloc = [];
else
    E.chromloc = Ein.chromloc;
end
if(~isfield(Ein,'annotation'))
    E.annotation = struct;
else
    E.annotation = Ein.annotation;
end
if(~isfield(Ein,'genegroups'))
    E.genegroups = cell(0);
else
    E.genegroups = Ein.genegroups;
end


% mache Länge von hybrname richtig
if(isempty(E.hybrname))
    for i=1:length(E.arrays)
        E.hybrname{i} = get(E.arrays{i},'name');
    end
end
% mache Länge von IDs richtig
if(isempty(E.IDs))
    for i=1:get(E.arrays{1},'ngene')
        E.IDs{i} = ['ID',num2str(i)];
    end
end

E = sortfields(E);
