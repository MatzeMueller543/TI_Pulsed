clearvars;clc;close all
RoundRoubin_ReferencePulsed

%% Aufgabe 1 Teil1: Fehlerwerte rausrechnen (9999 in Lidar_N.RWS bzw Lidar_S.RWS)

Fehlerwerte = Lidar_N.RWS==9999; %Logical Array 
%Test1=Lidar_N.RWS(Fehlerwerte)     Enthält alle 9999er Fehler (True = 1 - Werte)
%Test2=Lidar_N.RWS(~Fehlerwerte)    Enthält alle anderen Punkte außer 9999er Fehler (False = 0 - Werte)
Lidar_N.RWS(Fehlerwerte)=interp1(Lidar_N.t(~Fehlerwerte),Lidar_N.RWS(~Fehlerwerte),Lidar_N.t(Fehlerwerte)); %Korrektur der Fehlerwerte

Fehlerwerte = Lidar_S.RWS==9999;
Lidar_S.RWS(Fehlerwerte)=interp1(Lidar_S.t(~Fehlerwerte),Lidar_S.RWS(~Fehlerwerte),Lidar_S.t(Fehlerwerte)); %Korrektur der Fehlerwerte

%Zeitspanne zuweisen und in Variable schreiben next dann RWS und Mast bezogen auf Zeitvektor über
%Zeit plotten 
Zeitwerte_Lidar=datetime(2020,09,03,19,0,0:4:3600)';
Zeit_ref19bis20=Mast_N.TIMESTAMP(1:3601); %von Start Aufnahme bis Zeile 3601 -> 19 Uhr bis 20 Uhr am 

subplot(2,1,1);
plot(Zeit_ref19bis20,Reference.LOS_N(1:3601)) %Plot für Nord
xticks(Mast_N.TIMESTAMP(1:301:3601))
ylabel('LOS in m/s')
title('Windgeschwigkeit Nord')
hold on
plot(Zeitwerte_Lidar,Lidar_N.RWS(1:901))
legend('Reference_N', 'Lidar_N')
hold off

subplot(2,1,2)
plot(Zeit_ref19bis20,Reference.LOS_S(1:3601)) %Plot für Süd
xticks(Mast_S.TIMESTAMP(1:300:3601))
ylabel('LOS in m/s') 
title('Windgeschwigkeit Süd')
hold on
plot(Zeitwerte_Lidar,Lidar_S.RWS(1:901))
legend('Reference_S', 'Lidar_S')
hold off

%% Aufgabe 2 einfacher Vergleich pulsed (10min)
%%% get 10 min mean for Lidar

Timestamps_korr=datetime(2020,09,03,19,0,0:4:86400)';

t1=Timestamps_korr(1); %schleife geht hier los
t2=Timestamps_korr(end); %schleife endet hier

t_vec = datenum(t1:minutes(10):t2); % create a ideal time vector

n_10min = length(t_vec)-1; %143x 10 min am Tag

for i_10min= 1:n_10min
    Considered_N    	= Lidar_N.t>=t_vec(i_10min) & Lidar_N.t<t_vec(i_10min+1);
    Considered_S      	= Lidar_S.t>=t_vec(i_10min) & Lidar_S.t<t_vec(i_10min+1);  
    
    Lidar_10min.LOS_N_mean(i_10min) = nanmean(Lidar_N.RWS(Considered_N));
    Lidar_10min.LOS_S_mean(i_10min) = nanmean(Lidar_S.RWS(Considered_S));
    Lidar_10min.LOS_N_std(i_10min)  = nanstd (Lidar_N.RWS(Considered_N));
    Lidar_10min.LOS_S_std(i_10min)  = nanstd (Lidar_S.RWS(Considered_S));
end

Lidar_10min.LOS_TI_N = Lidar_10min.LOS_N_std./Lidar_10min.LOS_N_mean; %TI Berechnen Nord
Lidar_10min.LOS_TI_S = Lidar_10min.LOS_S_std./Lidar_10min.LOS_S_mean; %TI Berechnen Süd


p_N = polyfit(Reference_10min.LOS_TI_N,Lidar_10min.LOS_TI_N,1); %Regression (p1=m=Steigung p2=b=Y-Abschnitt)
p_S= polyfit(Reference_10min.LOS_TI_S,Lidar_10min.LOS_TI_S,1);


y_N=p_N(1)*Reference_10min.LOS_TI_N+p_N(2); %Geradengleichung Regressionsgerade
y_S=p_S(1)*Reference_10min.LOS_TI_S+p_S(2);

r_N = corrcoef(Reference_10min.LOS_TI_N,Lidar_10min.LOS_TI_N); %R Wert (R^2 = Bestimmtheitsmaß)
r_S = corrcoef(Reference_10min.LOS_TI_S,Lidar_10min.LOS_TI_S);

r_srt_N = r_N(1,2)^2 ;%Bestimmtsheitsmaß
r_srt_N_string = ['R^2 = ' , num2str(r_srt_N)]; %zu String für Plot

r_srt_S = r_S(1,2)^2;
r_srt_S_string = ['R^2 = ' , num2str(r_srt_S)]; %zu String für Plot
 
subplot(2,1,1);
plot(Reference_10min.LOS_TI_N,Lidar_10min.LOS_TI_N,'b.')
text(0.5,0.24,r_srt_N_string);
title('Lidar TI_N über Reference TI_N')
grid on; grid minor
xlabel('TI Reference')
ylabel('TI Lidar')
hold on
plot(Reference_10min.LOS_TI_N,y_N) %plotte Regressionsgerade
hold off

subplot(2,1,2);
plot(Reference_10min.LOS_TI_S,Lidar_10min.LOS_TI_S,'b.')
text(0.35,0.13,r_srt_S_string);
title('Lidar TI_S über Reference TI_S')
grid on; grid minor
xlabel('TI Reference')
ylabel('TI Lidar')
hold on
plot(Reference_10min.LOS_TI_S,y_S) %plotte Regressionsgerade
hold off














