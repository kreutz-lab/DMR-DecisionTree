% Konstruktor der Klasse microarray
% M = microarray(dat, name)

function M = microarray(dat, name);
if(~exist('name','var') | isempty(name))
	name = '';
end

if(~isstruct(dat))
	M.name  = name;
	M.data  = dat;
	
	M.discreteFactors   = [];
	M.continuousFactors = [];
	
	M = class(M,'microarray');
else
    M = class(dat,'microarray');
end
