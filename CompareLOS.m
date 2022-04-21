function CompareLOS(Reference_10min,Lidar_10min)

range_WS = [0 10];
range_TI = [0 2.5];

figure('Name','LOS mean and std. regression')

subplot(2,2,1)
hold on;box on;grid on
plot(Reference_10min.LOS_N_mean,Lidar_10min.LOS_N_mean,'.');
plot(range_WS,range_WS);
title('Mean LOS, North')
xlabel('Reference [m/s]')
ylabel('Lidar [m/s]')
axis equal
xlim(range_WS)
ylim(range_WS)
p = polyfit(Reference_10min.LOS_N_mean,Lidar_10min.LOS_N_mean,1);
text(0.1*range_WS(2),0.9*range_WS(2),['y=',num2str(p(2),'%4.2f'),'+',num2str(p(1),'%4.2f'),' x'])

subplot(2,2,2)
hold on;box on;grid on
plot(Reference_10min.LOS_S_mean,Lidar_10min.LOS_S_mean,'.')
plot(range_WS,range_WS);
title('Mean LOS, South')
xlabel('Reference [m/s]')
ylabel('Lidar [m/s]')
axis equal
xlim(range_WS)
ylim(range_WS)
p = polyfit(Reference_10min.LOS_S_mean,Lidar_10min.LOS_S_mean,1);
text(0.1*range_WS(2),0.9*range_WS(2),['y=',num2str(p(2),'%4.2f'),'+',num2str(p(1),'%4.2f'),' x'])


subplot(2,2,3)
hold on;box on;grid on
plot(Reference_10min.LOS_N_std,Lidar_10min.LOS_N_std,'.')
plot(range_TI,range_TI);
title('STD LOS, North')
xlabel('Reference [m/s]')
ylabel('Lidar [m/s]')
axis equal
xlim(range_TI)
ylim(range_TI)
p = polyfit(Reference_10min.LOS_N_std,Lidar_10min.LOS_N_std,1);
text(0.1*range_TI(2),0.9*range_TI(2),['y=',num2str(p(2),'%4.2f'),'+',num2str(p(1),'%4.2f'),' x'])

subplot(2,2,4)
hold on;box on;grid on
plot(Reference_10min.LOS_S_std,Lidar_10min.LOS_S_std,'.')
plot(range_TI,range_TI);
title('STD LOS, South')
xlabel('Reference [m/s]')
ylabel('Lidar [m/s]')
axis equal
xlim(range_TI)
ylim(range_TI)
p = polyfit(Reference_10min.LOS_S_std,Lidar_10min.LOS_S_std,1);
text(0.1*range_TI(2),0.9*range_TI(2),['y=',num2str(p(2),'%4.2f'),'+',num2str(p(1),'%4.2f'),' x'])


end