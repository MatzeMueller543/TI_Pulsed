function Coordinate =  CalculateCoordinate
% function to calculate all necessary coordinates into inertial coordinate
% system: Origin is location of lidar system, x to east, y to north)

% coordinates in UTM32 (easting, northing and height above ground)
UTM.Mast_N                  = [520371.82,6054567.54,30];
UTM.Mast_S                  = [520433.26,6054463.81,30];
UTM.Lidar                   = [520555.96,6054607.73,30]; 

% reduce to Inertial coordinate system (origin at lidar system location)
Coordinate.Mast_N           = UTM.Mast_N-UTM.Lidar;   
Coordinate.Mast_S           = UTM.Mast_S-UTM.Lidar;   
Coordinate.Lidar            = [0,0,0];   

% get yaw angles of met masts
[YawMast_N,~,~]             = cart2sph(Coordinate.Mast_N(1),Coordinate.Mast_N(2),Coordinate.Mast_N(3));
[YawMast_S,~,~]             = cart2sph(Coordinate.Mast_S(1),Coordinate.Mast_S(2),Coordinate.Mast_S(3));

% scan configuration
VerticalHalfOpeningAngle    = deg2rad(5);
HorizontalHalfOpeningAngle  = deg2rad(15);
MeasurementDistance         = 178;

% coordinates of the 2 lower beams in the Lidar coordinate system
x_F_L                       = MeasurementDistance;
y_F_L                       = x_F_L.*[-1  1]*tan(HorizontalHalfOpeningAngle);
z_F_L                       = x_F_L.*[-1 -1]*tan(VerticalHalfOpeningAngle);

% transformation into inertial coordinate system
YawLidar                    = (YawMast_N+YawMast_S)/2+deg2rad(1.6);
PitchLidar                  = -deg2rad(5);
RollLidar                   = 0;
[x_F_I,y_F_I,~]             = L2I(x_F_L,y_F_L,z_F_L,YawLidar,PitchLidar,RollLidar);

% store into Cordinate struct
Coordinate.Focus_N          = [x_F_I(1),y_F_I(1),0];
Coordinate.Focus_S          = [x_F_I(2),y_F_I(2),0]; 

end


function [x_I,y_I,z_I] = L2I(x_L,y_L,z_L,Yaw,Pitch,Roll)

% Yaw is a rotation around z-axis    
T_Yaw 	= [ cos(Yaw)   -sin(Yaw)	0;
            sin(Yaw) 	cos(Yaw)    0;
            0           0           1];

% Pitch is a rotation around y-axis
T_Pitch = [	cos(Pitch)	0           sin(Pitch);
            0         	1        	0;
           -sin(Pitch)	0        	cos(Pitch)];

% Roll is a rotation around x-axis
T_Roll  = [	1           0       	0;
            0       	cos(Roll)  -sin(Roll);
            0          	sin(Roll)	cos(Roll)];

% Rotation Order is yaw-pitch-roll        
T       = T_Yaw*T_Pitch*T_Roll;

% Transformation
x_I     = T(1,1)*x_L + T(1,2)*y_L + T(1,3)*z_L;
y_I     = T(2,1)*x_L + T(2,2)*y_L + T(2,3)*z_L;
z_I     = T(3,1)*x_L + T(3,2)*y_L + T(3,3)*z_L;

end