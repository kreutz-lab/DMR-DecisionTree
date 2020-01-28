% out = imagesc_nan(varargin)
% 
%   imagesc, alle NaN werden auf -inf gesetzt und weiss geplottet.

function out = imagesc_nan(dat,cmap,varargin)
if nargin<2
    cmap = redgreencmap;
end
imagesc(dat,varargin{:});
c = colormap(cmap);
c = c(end:-1:1,:);
c = [ones(1,3);c];
dat(isnan(dat)) = -.05;
% varargin{1} = varargin{1}+1;
% varargin{1}(varargin{1}==0) = ;
out = imagesc(dat,varargin{:});
colormap(c)


