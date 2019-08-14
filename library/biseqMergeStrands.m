% dat = biseqMergeStrands(dat)
% 
%   Merging of reads for forward and backward position
% 
% Data from forward and backward strands are merged.
% Methylation is symetric, i.e. should occur on both strands. The negative
% strand position is usually 1 larger than the respective position on the
% postive strand. 
% 
%   If only a backward positon is available, then the number of reads at
%   the respective forward positon was zero.


function dat = biseqMergeStrands(dat)

if length(unique(dat.isF)) ~=2
    unique(dat.isF)
    error('length(unique(dat.isF))~=2: Either the strands were already merged or the strand info is not available.')
end

datIn = dat;

chrlev = unique(dat.chr);
raus = [];
for i=1:length(chrlev)
    indchr = find(ismember(dat.chr,chrlev{i}));
    
    isF = logical(dat.isF(indchr)); 
    indIsF = find(isF);
    pos = dat.pos(indchr);
    posF = pos(isF);
    indNotF = find(~isF);
    posB = pos(indNotF);
    
    [~,indF,indB] = intersect(posF,posB-1);
    if ~isempty(indF)
        fprintf('%i (out of %i) positions from forward and backward strands are merged for chromosome %s\n',length(indF),length(dat.pos),chrlev{i})
    end
    dat.nread(indchr(indIsF(indF))) = dat.nread(indchr(indIsF(indF))) + dat.nread(indchr(indNotF(indB))); % sum up the reads        
    
%     if length(indF)>100
%         plot(dat.nread(indchr(indF))+randn(size(indF))*.3,dat.nread(indchr(indB))+randn(size(indB))*.3,'.')
%         title(['Agreement of data from forward and backward stands for chr',str2label(chrlev{i})]) 
%         xlabel('nread+randn forward strand')
%         ylabel('nread+randn backward strand');
%         abplot(1,0);
% %         print(gcf,['biseqMergesSymmetry_',chrlev{i}],'-dpng')
%     end
    
    raus = [raus;indchr(indNotF(indB))]; 
    
    
    % now handle pos with only reads on the backward strand:
    [~,indB] = setdiff(posB-1,posF);
    dat.pos(indchr(indNotF(indB))) = dat.pos(indchr(indNotF(indB)))-1; % sum up the reads        
    for j=1:length(indB)
        dat.strand{indchr(indNotF(indB(j)))} = '0';
    end
end

drin = 1:length(dat.pos);
drin(raus) = [];
% fprintf('A total of %i reads were merged out of %i positions.\n',length(raus),length(dat.pos));
dat = FilterStruct(dat,drin,length(dat.pos),1);


%% do some checks:
% if sum(datIn.nread) ~= sum(dat.nread) % total number of reads should be constant
%     fprintf('%i reads from forward and backward strands are merged.',length(raus));
%     error('sum(datIn.nread) ~= sum(dat.nread)');
% end
if(sum(datIn.isF)>sum(dat.isF)) % % no forward positions should disappear
    error('sum(datIn.isF)>sum(dat.isF)')
end
if(sum(~datIn.isF) - sum(~dat.isF) ~= length(raus)) % the number of backward positions which disappear should conicide with the numer of removed positions
    error('sum(~datIn.isF) - sum(~dat.isF) ~= length(raus)')
end



