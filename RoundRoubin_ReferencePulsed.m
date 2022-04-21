% Calculate Reference TI
% Project: IEA Wind Task 32
% Round Robin on turbulence estimates from nacelle mounted lidar systems
% by Feng Guo and David Schlipf @ Flensburg University of Applied Sciences
% v4: 11-Apr-2022: add new data
% v3: 10-Apr-2022: store data to speed up 2nd run 
% v2: 14-Mar-2022: modular setup
% v1: 30-Nov-2021: initial version

clearvars;clc;close all
%addpath('functions')

%% Calculate and plot coordinates
Coordinate      =  CalculateCoordinate;
PlotCoordinate(Coordinate)

%% Load Data and compare Sonic Cup and Vane
if isfile('Data.mat') % datenum takes a while, so we better store the data
    load('Data.mat','Mast_N','Mast_S','Lidar_N','Lidar_S');
else 
    Mast_N          = readtable('TMMN_20200903_20200904_1Hz_new.csv'); % lines 14704-14707 with text removed 
    Mast_S          = readtable('TMMS_20200903_20200904_1Hz_new.csv'); % lines 14704-14707 with text removed 
    Lidar_N        	= readtable('Lidar_20200903-20200904_1Hz_LOS3_178m.csv');
    Lidar_S       	= readtable('Lidar_20200903-20200904_1Hz_LOS2_178m.csv');
    % remove double entries 
    [~,IdxUnique,~] = unique(Mast_S.TIMESTAMP);
    Mast_S          = Mast_S(IdxUnique,:);  
    % add numeric time
    Mast_N.t      	= datenum(Mast_N.TIMESTAMP);
    Mast_S.t     	= datenum(Mast_S.TIMESTAMP);
    Lidar_N.t     	= datenum(Lidar_N.Timestamp,'yyyy-mm-ddTHH:MM:SS.FFF');
    Lidar_S.t      	= datenum(Lidar_S.Timestamp,'yyyy-mm-ddTHH:MM:SS.FFF');
    save('Data.mat','Mast_N','Mast_S','Lidar_N','Lidar_S');
end

% Remove manually some outliers for cup and interpolate
GoodData            = [1:82155,82210:82260,82310:86400];% identified using cursor
Mast_N.WS1          = interp1(Mast_N.t(GoodData),Mast_N.WS1(GoodData),Mast_N.t);

CompareSonicToCupAndVane(Mast_S,Mast_N)

%% Calculate reference U (wind in x: W->E), V (wind in y: S->N)
Yaw_S              	= 270-Mast_S.USA_WD;
Yaw_N             	= 270-Mast_N.USA_WD;
Reference.U_S      	= Mast_S.WS1.*cosd(Yaw_S);
Reference.V_S     	= Mast_S.WS1.*sind(Yaw_S);
Reference.U_N    	= Mast_N.WS1.*cosd(Yaw_N);
Reference.V_N     	= Mast_N.WS1.*sind(Yaw_N);
Reference.t         = Mast_N.t; % should be the same for N and S

%% Calculate reference LOS
% unit vectors
LidarVector_N    	= [Coordinate.Focus_N(1) Coordinate.Focus_N(2) Coordinate.Focus_N(3)];
LidarVector_S    	= [Coordinate.Focus_S(1) Coordinate.Focus_S(2) Coordinate.Focus_S(3)];
LidarUnitVector_N 	= LidarVector_N/norm(LidarVector_N);
LidarUnitVector_S  	= LidarVector_S/norm(LidarVector_S);
% measurement equation: ignoring w, convention LOS is positive if towards system
Reference.LOS_S    	= -(LidarUnitVector_S(1)*Reference.U_S + LidarUnitVector_S(2)*Reference.V_S);
Reference.LOS_N    	= -(LidarUnitVector_N(1)*Reference.U_N + LidarUnitVector_N(2)*Reference.V_N);

%% get 10 min mean for reference 
Tstart          = '2020-09-03 19:00:00';
Tend            = '2020-09-04 19:00:00';
Reference_10min = Calculate10minStastics_Reference(Reference,Tstart,Tend);

%% get 10 min reference TI
Reference_10min.LOS_TI_N = Reference_10min.LOS_N_std./Reference_10min.LOS_N_mean;
Reference_10min.LOS_TI_S = Reference_10min.LOS_S_std./Reference_10min.LOS_S_mean;
WriteReferenceFile(Reference_10min,'ReferencePulsed.csv')
Reference.U_S;
