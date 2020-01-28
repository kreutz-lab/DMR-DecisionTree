% This function extracts parameters used as predictors from simulation arg
% struct.
% The selection was made based on the assignmen in wgbss_optimize_p2arg.m
% 
%   donorm  controls whether the values val are standardized


function [val, names] = args2predictor(args,donorm)
if ~exist('donorm','var') || isempty(donorm)
    donorm = true;
end

val = NaN(length(args.arg),14);
names = cell(1,14);
names{1} = 'P(succ in meth)';
names{2} = 'P(succ in nonmeth)';
names{3} = 'Err(meth)';
names{4} = 'Err(nonmeth)';
names{5} = 'Mean(meth)';
names{6} = 'Mean(nonmeth)';
names{7} = 'Dense(Island)';
names{8} = 'Dense(Desert)';
names{9} = 'Pdecay';
names{10} = 'StayInIsland';
names{11} = 'LeaveDesert';
names{12} = 'LeaveMeth';
names{13} = 'LeaveNonmeth';
names{14} = 'invertT';

for i=1:length(args.arg)
    arg = args.arg{i};

    val(i,1) = arg.probability_of_success_in_methylated_region;
    val(i,2) = arg.probability_of_success_in_non_methylated_region;
    val(i,3) = arg.error_rate_in_methylated_region;
    val(i,4) = arg.error_rate_in_non_methylated_region;
    val(i,5) = arg.mean_number_of_reads_in_methylated_region;
    val(i,6) = arg.mean_number_of_reads_in_non_methylated_region;
    
    val(i,7) = arg.cpg_matrix{1}(2);
    val(i,8) = arg.cpg_matrix{1}(1);
    val(i,9) = arg.dist_value;
    val(i,10) = arg.transCpgLoc(3);
    val(i,11) = arg.transCpgLoc(1);
    
    if isfield(arg,'invertT') && arg.invertT  % invertT=true is proper interpretation/implementation
        val(i,7) = arg.cpg_matrix{1}(1); % invertT seems anticorrelated to these two parameters
        val(i,8) = arg.cpg_matrix{1}(2);
        
        val(i,10) = arg.transCpgLoc(1); % invertT changes the meaning
        val(i,11) = arg.transCpgLoc(3);
    else
        arg.invertT = false;
    end
    val(i,12) = arg.transPi{1}(1);
    val(i,13) = arg.transPi{1}(3);
    val(i,14) = arg.invertT;
end

val = real(val); % eine Zahl ist imaginär, keine Ahnung wieso

% if donorm
%     for i=1:size(val,2)
%         val(:,i) = val(:,i)-mean(val(:,i));
%         val(:,i) = val(:,i)./std(val(:,i));
%     end
% end
