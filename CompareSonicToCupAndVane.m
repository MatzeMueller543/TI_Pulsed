function CompareSonicToCupAndVane(Mast_S,Mast_N)

figure('Name','Comparison over time')

subplot(2,2,1)
hold on;box on;grid on
plot(sqrt(Mast_S.USA_U.^2+Mast_S.USA_V.^2))
plot(Mast_S.WS1)
title('South Mast: horizontal wind speed')
ylabel('[m/s]')
legend('Sonic','Cup')

subplot(2,2,2)
hold on;box on;grid on
plot(sqrt(Mast_N.USA_U.^2+Mast_N.USA_V.^2))
plot(Mast_N.WS1)
title('North Mast: horizontal wind speed')
ylabel('[m/s]')
legend('Sonic','Cup')

subplot(2,2,3)
hold on;box on;grid on
plot(Mast_S.USA_WD)
plot(Mast_S.WD1)
title('South Mast: wind direction')
ylabel('[deg]')
xlabel('data point [-]')
legend('Sonic','Vane')

subplot(2,2,4)
hold on;box on;grid on
plot(Mast_N.USA_WD)
plot(Mast_N.WD1)
title('North Mast: wind direction')
ylabel('[deg]')
xlabel('data point [-]')
legend('Sonic','Vane')


figure('Name','Regression')

subplot(2,2,1)
hold on;box on;grid on
plot(sqrt(Mast_S.USA_U.^2+Mast_S.USA_V.^2),Mast_S.WS1,'.')
title('South Mast: horizontal wind speed')
xlabel('Sonic [m/s]')
ylabel('Cup [m/s]')

subplot(2,2,2)
hold on;box on;grid on
plot(sqrt(Mast_N.USA_U.^2+Mast_N.USA_V.^2),Mast_N.WS1,'.')
title('North Mast: horizontal wind speed')
xlabel('Sonic [m/s]')
ylabel('Cup [m/s]')

subplot(2,2,3)
hold on;box on;grid on
plot(Mast_S.USA_WD,Mast_S.WD1,'.')
title('South Mast: wind direction')
xlabel('Sonic [deg]')
ylabel('Cup [deg]')

subplot(2,2,4)
hold on;box on;grid on
plot(Mast_N.USA_WD,Mast_N.WD1,'.')
title('North Mast: wind direction')
xlabel('Sonic [deg]')
ylabel('Cup [deg]')

end
