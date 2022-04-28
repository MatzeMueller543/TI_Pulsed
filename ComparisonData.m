function ComparisonData(Reference_10min,Lidar_10min)

p_R_TI = polyfit(Reference_10min.LOS_TI_N,Reference_10min.LOS_TI_S,1);
p_R_mean = polyfit(Reference_10min.LOS_N_mean,Reference_10min.LOS_S_mean,1);
p_R_std = polyfit(Reference_10min.LOS_N_std,Reference_10min.LOS_S_std,1);
p_L_TI = polyfit(Lidar_10min.LOS_TI_N,Lidar_10min.LOS_TI_S,1);
p_L_mean = polyfit(Lidar_10min.LOS_N_mean,Lidar_10min.LOS_S_mean,1);
p_L_std = polyfit(Lidar_10min.LOS_N_std,Lidar_10min.LOS_S_std,1);

x_max = [0, 0.5];
x_max_2 = [0, 10];

y_R_TI = polyval(p_R_TI,x_max);
y_R_mean = polyval(p_R_mean,x_max_2);
y_R_std = polyval(p_R_std,x_max_2);
y_L_TI = polyval(p_L_TI,x_max);
y_L_mean = polyval(p_L_mean,x_max_2);
y_L_std = polyval(p_L_std,x_max_2);

r_R_TI = corrcoef(Reference_10min.LOS_TI_N,Reference_10min.LOS_TI_S); 
r_R_mean = corrcoef(Reference_10min.LOS_N_mean,Reference_10min.LOS_S_mean); 
r_R_std = corrcoef(Reference_10min.LOS_N_std,Reference_10min.LOS_S_std); 
r_L_TI = corrcoef(Lidar_10min.LOS_TI_N,Lidar_10min.LOS_TI_S);
r_L_mean = corrcoef(Lidar_10min.LOS_N_mean,Lidar_10min.LOS_S_mean);
r_L_std = corrcoef(Lidar_10min.LOS_N_std,Lidar_10min.LOS_S_std);

r_sq_R_TI = r_R_TI(1,2)^2;
r_sq_R_str_TI = ['R^2 = ' , num2str(r_sq_R_TI)];
r_sq_R_mean = r_R_mean(1,2)^2;
r_sq_R_str_mean = ['R^2 = ' , num2str(r_sq_R_mean)];
r_sq_R_std = r_R_std(1,2)^2;
r_sq_R_str_std = ['R^2 = ' , num2str(r_sq_R_std)];
r_sq_L_TI = r_L_TI(1,2)^2;
r_sq_L_str_TI = ['R^2 = ' , num2str(r_sq_L_TI)];
r_sq_L_mean = r_L_mean(1,2)^2;
r_sq_L_str_mean = ['R^2 = ' , num2str(r_sq_L_mean)];
r_sq_L_std = r_L_std(1,2)^2;
r_sq_L_str_std = ['R^2 = ' , num2str(r_sq_L_std)];

x_r_pos = 0.5*x_max(2);
y_r_pos = 0.2*x_max(2);

figure('name','comparison Referencedata');
subplot(2,2,1);
hold on; grid on;box on;
plot(Reference_10min.LOS_N_mean,Reference_10min.LOS_S_mean,'.')
plot(Reference_10min.LOS_N_std,Reference_10min.LOS_S_std,'.')
plot(x_max_2,y_R_mean)
plot(x_max_2,y_R_std)
title('Reference mean and std N to S')
xlabel('mean and std Reference_N')
ylabel('mean and std Reference_S')
axis equal
xlim(x_max_2)
ylim(x_max_2)

subplot(2,2,2);
hold on; grid on;box on;
plot(Lidar_10min.LOS_N_mean,Lidar_10min.LOS_S_mean,'.')
plot(Lidar_10min.LOS_N_std,Lidar_10min.LOS_S_std,'.')
plot(x_max_2,y_L_mean)
plot(x_max_2,y_L_std)
title('Lidar mean and std N to S')
xlabel('mean and std Lidar_N')
ylabel('mean and std Lidar_S')
axis equal
xlim(x_max_2)
ylim(x_max_2)

subplot(2,2,3);
hold on; grid on;box on;
plot(Reference_10min.LOS_TI_N,Reference_10min.LOS_TI_S,'.')
plot(x_max,y_R_TI)
text(x_r_pos,y_r_pos,r_sq_R_str_TI); 
title('Reference TI_S über Reference TI_N')
xlabel('TI Reference_N')
ylabel('TI Reference_S')
axis equal
xlim(x_max)
ylim(x_max)

subplot(2,2,4);
hold on; grid on;box on;
plot(Lidar_10min.LOS_TI_N,Lidar_10min.LOS_TI_S,'.')
plot(x_max,y_L_TI)
text(x_r_pos,y_r_pos,r_sq_L_str_TI); 
title('Lidar TI_S über Lidar TI_N')
xlabel('TI Lidar_N')
ylabel('TI Lidar_S')
axis equal
xlim(x_max)
ylim(x_max)

end

