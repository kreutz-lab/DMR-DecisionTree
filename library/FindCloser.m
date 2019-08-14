%   FindCloser(bol1,bol2)
%   FindCloser(bol1,bol2,x)
% 
%   Findet für alle indices ind ob bol1 oder bol2 näher ist
% 
%   x  positions of bol1 and bol2
% 
% ind   array of indices
% bol1  logical array indicating positions of type1
% bol2  logical array indicating positions of type2
% 
%   indicator   1 if bol1==true näher
%               2 if bol2==true näher
%               1.5 if both equidistant

function indicator = FindCloser(bol1,bol2,x)
if ~exist('x','var') || isempty(x)
    dx = ones(size(bol1));
else   
    dx = diff(x);
end

if size(bol1,1)>1 
    bol1 = bol1';
end
if size(bol2,1)>1 
    bol2 = bol2';
end

if sum(abs(size(bol1)-size(bol2)))~=0
    error('FindCloser.m: size(bol1) unequal to size(bol2)')
end

valleft1 = double(bol1);
for i=2:length(bol1)
    if ~bol1(i)
        valleft1(i) = valleft1(i-1)*(.999^dx(i-1));
    end
end

valleft2 = double(bol2);
for i=2:length(bol2)
    if ~bol2(i)
        valleft2(i) = valleft2(i-1)*(.999^dx(i-1));
    end
end

valright1 = double(bol1);
for i=(length(bol1)-1):-1:1
    if ~bol1(i)
        valright1(i) = valright1(i+1)*(.999^dx(i));
    end
end

valright2 = double(bol2);
for i=(length(bol2)-1):-1:1
    if ~bol2(i)
        valright2(i) = valright2(i+1)*(.999^dx(i));
    end
end

formax = [valleft1;valright1;valleft2;valright2];
[~,ind] = max(formax);
[~,ind2] = max(flip(formax,1));
ind2 = size(formax,1)+1-ind2;

indicator = 0.5*ceil(ind/2) + 0.5*ceil(ind2/2); % average bei gleichem Abstand
