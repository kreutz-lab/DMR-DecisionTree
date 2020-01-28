%PRINTALLTOPng macht aus mehreren Figures in ein Png-file
%
%PRINTALLTOPng(fignr,dateiname) Schreibt figure 1:fignr
%	in dateiname.jpg
%
%PRINTALLTOPng(fignr,dateiname,opt) 
%   opt=1   Figures schlieﬂen
%
function PrintToPng(fignr,dateiname,opt)
if(~exist('fignr','var'))
	fignr=[];
end
disp('Png wird geschrieben');
set(fignr,'InvertHardCopy','off','Color',ones(1,3)); % fur Hintergrund der subplots/axes

try
    print(fignr,'-dpng',[dateiname,'.png']);
catch
    fignr
    dateiname
    disp(lasterr)
end

if(exist('opt','var'))
    if (opt==1)
        close all;
    end
end