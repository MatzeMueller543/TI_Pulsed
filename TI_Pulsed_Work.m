clearvars;clc;close all
RoundRoubin_ReferencePulsed

%% Task 1 Part1: clean data (9999 in Lidar_N.RWS bzw Lidar_S.RWS)

Fehlerwerte = Lidar_N.RWS==9999; %Logical Array 
Lidar_N.RWS(Fehlerwerte)=interp1(Lidar_N.t(~Fehlerwerte),Lidar_N.RWS(~Fehlerwerte),Lidar_N.t(Fehlerwerte)); 

Fehlerwerte = Lidar_S.RWS==9999;
Lidar_S.RWS(Fehlerwerte)=interp1(Lidar_S.t(~Fehlerwerte),Lidar_S.RWS(~Fehlerwerte),Lidar_S.t(Fehlerwerte)); 

MyXlim=[datenum('2020-09-03 19:00:00') datenum('2020-09-03 20:00:00')];

figure('name','time comparison windspeed')
subplot(2,1,1);
hold on; box on; grid on;
plot(Reference.t,Reference.LOS_N)
plot(Lidar_N.t,Lidar_N.RWS)
xlim(MyXlim)
datetick('x','keeplimits')
ylabel('LOS in m/s') 
title('Windgeschwigkeit Nord')
legend('Reference_N', 'Lidar_N')

subplot(2,1,2);
hold on; box on; grid on;
plot(Reference.t,Reference.LOS_S)
plot(Lidar_S.t,Lidar_S.RWS)
xlim(MyXlim)
datetick('x','keeplimits')
ylabel('LOS in m/s') 
title('Windgeschwigkeit Süd')
legend('Reference_S', 'Lidar_S')

%% Task 2 comparison pulsed (simple 10min)

Lidar_10min = Calculate10minStastics_Lidar(Lidar_N,Lidar_S,Tstart,Tend);


Lidar_10min.LOS_TI_N = Lidar_10min.LOS_N_std./Lidar_10min.LOS_N_mean; %TI North
Lidar_10min.LOS_TI_S = Lidar_10min.LOS_S_std./Lidar_10min.LOS_S_mean; %TI South


p_N = polyfit(Reference_10min.LOS_TI_N,Lidar_10min.LOS_TI_N,1); %Regression
p_S = polyfit(Reference_10min.LOS_TI_S,Lidar_10min.LOS_TI_S,1);

x_N = [0 , 0.7];
x_S = x_N;
y_N = polyval(p_N, x_N); 
y_S = polyval(p_S, x_S);

r_N = corrcoef(Reference_10min.LOS_TI_N,Lidar_10min.LOS_TI_N); %R = corealcoefficient ()
r_S = corrcoef(Reference_10min.LOS_TI_S,Lidar_10min.LOS_TI_S);

r_srt_N = r_N(1,2)^2 ;%R^2 = determination coefficient
r_srt_N_string = ['R^2 = ' , num2str(r_srt_N)]; %zu String für Plot

r_srt_S = r_S(1,2)^2;
r_srt_S_string = ['R^2 = ' , num2str(r_srt_S)]; %zu String für Plot

x_r = 0.5*x_N(2);
y_r = 0.2*x_N(2);

figure('name','comparison TI');
subplot(1,2,1);
hold on; grid on; box on;
plot(x_N,y_N) 
plot(Reference_10min.LOS_TI_N,Lidar_10min.LOS_TI_N,'b.')
text(x_r,y_r,r_srt_N_string); %text außerhalb
title('Lidar TI_N über Reference TI_N')
xlabel('TI Reference')
ylabel('TI Lidar')
axis equal
xlim(x_N)
ylim(x_N)

subplot(1,2,2);
hold on; grid on; box on;
plot(x_S,y_S) 
plot(Reference_10min.LOS_TI_S,Lidar_10min.LOS_TI_S,'b.')
text(x_r,y_r,r_srt_S_string);
title('Lidar TI_S über Reference TI_S')
xlabel('TI Reference')
ylabel('TI Lidar')
axis equal
xlim(x_S)
ylim(x_S)

%% Comparison Reference and Lidar data
ComparisonData(Reference_10min,Lidar_10min) 



