% Renaming the data context (sample names) according to Noe's suggestiong
% (3.6.19)

function in = ReplaceSampleNames(in)

in = strrep(in,'P.ab-','Picab-');
in = strrep(in,'PP-','Phypa-');
in = strrep(in,'A.th-','Arath-');
in = strrep(in,'A.th ','Arath ');
in = strrep(in,'A.th_','Arath_');
in = strrep(in,'Ae-','Aetar-');

in = strrep(in,'CpG','CG');

in = strrep(in,'moabs','MOABS');
in = strrep(in,'dmrcate','DMRcate');
in = strrep(in,'defiant','Defiant');
in = strrep(in,'bsmooth','BSmooth');
in = strrep(in,'methylsig','MethylSig');
in = strrep(in,'methylkit','MethylKit');
in = strrep(in,'methylscore','MethylScore');



