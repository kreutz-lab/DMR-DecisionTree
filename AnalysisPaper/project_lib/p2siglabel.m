% Ermittelt zu p-Werten die *** label.
% 
%   thresholds: default [0.05,0.01,0.001]
%   d.h. * f�r p in [0.01,...,0.05]
%       ** f�r p in [0.001,...,0.01]
%      *** f�r p <0.001

function sl = p2siglabel(p,thresholds)

if (~exist('thresholds','var') || isempty(thresholds))
    thresholds = [0.05,0.01,0.001];
end

sl = cell(size(p));
for i1=1:size(p,1)
    for i2=1:size(p,2)
        sl{i1,i2} = '';
        for it = 1:length(thresholds)
            if(p(i1,i2)<thresholds(it))
                sl{i1,i2} = [sl{i1,i2},'*'];
            end
        end
        if(length(sl{i1,i2})==1)
            sl{i1,i2} = ' * ';
        end
    end
end

