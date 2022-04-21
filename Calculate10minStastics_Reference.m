function Reference_10min = Calculate10minStastics_Reference(Reference,Tstart,Tend) 
% function to calculate the 10min statistic for the Mast

% build time vector
t1                          = datetime(Tstart); 
t2                          = datetime(Tend); 
t_vec                       = datenum(t1:minutes(10):t2); % create a ideal time vector 

% loop over the data and calculate the 10min statistic
n_10min                     = length(t_vec)-1;
for i_10min= 1:n_10min
    Considered            = Reference.t>=t_vec(i_10min) & Reference.t<t_vec(i_10min+1);
    
    Reference_10min.LOS_N_mean(i_10min)   = nanmean(Reference.LOS_N(Considered));
    Reference_10min.LOS_S_mean(i_10min)   = nanmean(Reference.LOS_S(Considered));
    Reference_10min.LOS_N_std(i_10min)    = nanstd (Reference.LOS_N(Considered));
    Reference_10min.LOS_S_std(i_10min)    = nanstd (Reference.LOS_S(Considered));
    Reference_10min.U_N_mean(i_10min)     = nanmean(Reference.U_N  (Considered));
    Reference_10min.U_S_mean(i_10min)     = nanmean(Reference.U_S  (Considered));
    Reference_10min.V_N_mean(i_10min)     = nanmean(Reference.V_N  (Considered));
    Reference_10min.V_S_mean(i_10min)     = nanmean(Reference.V_S  (Considered));
end

% add time vector
Reference_10min.t           = t_vec(1:end-1);
Reference_10min.Timestamp   = datestr(t_vec(1:end-1),31);

end