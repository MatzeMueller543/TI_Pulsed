function Lidar_10min  = Calculate10minStastics_Lidar(Lidar_N,Lidar_S,Tstart,Tend)
% function to calculate the 10min statistic for the lidar

% build time vector
t1                    	= datetime(Tstart); 
t2                   	= datetime(Tend); 
t_vec                   = datenum(t1:minutes(10):t2); % create a ideal time vector 

% loop over the data and calculate the 10min statistic
n_10min                 = length(t_vec)-1;
for i_10min= 1:n_10min
    Considered_N    	= Lidar_N.t>=t_vec(i_10min) & Lidar_N.t<t_vec(i_10min+1);
    Considered_S      	= Lidar_S.t>=t_vec(i_10min) & Lidar_S.t<t_vec(i_10min+1);    
    Lidar_10min.LOS_N_mean(i_10min) = nanmean(Lidar_N.RWS(Considered_N));
    Lidar_10min.LOS_S_mean(i_10min) = nanmean(Lidar_S.RWS(Considered_S));
    Lidar_10min.LOS_N_std(i_10min)  = nanstd (Lidar_N.RWS(Considered_N));
    Lidar_10min.LOS_S_std(i_10min)  = nanstd (Lidar_S.RWS(Considered_S));
end

% add time vector
Lidar_10min.t           = t_vec(1:end-1);
Lidar_10min.Timestamp   = datestr(t_vec(1:end-1),31);

end