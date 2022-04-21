function PlotCoordinate(Coordinate)

figure('Name','Coordinates')
hold on;box on;grid on
plot3(Coordinate.Mast_N(1),Coordinate.Mast_N(2),Coordinate.Mast_N(3),'o')
plot3(Coordinate.Mast_S(1),Coordinate.Mast_S(2),Coordinate.Mast_S(3),'o')
plot3(Coordinate.Lidar(1), Coordinate.Lidar(2),Coordinate.Lidar(3),'^')
plot3([Coordinate.Lidar(1)  Coordinate.Focus_N(1)],[Coordinate.Lidar(2)  Coordinate.Focus_N(2)],[Coordinate.Lidar(3)  Coordinate.Focus_N(3)],'.-')
plot3([Coordinate.Lidar(1)  Coordinate.Focus_S(1)],[Coordinate.Lidar(2)  Coordinate.Focus_S(2)],[Coordinate.Lidar(3)  Coordinate.Focus_S(3)],'.-')
plot3([Coordinate.Mast_N(1) Coordinate.Mast_S(1) ],[Coordinate.Mast_N(2) Coordinate.Mast_S(2) ],[Coordinate.Mast_N(3) Coordinate.Mast_S(3) ],'--')
axis equal
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
legend('Mast N','Mast S','Lidar','Beam N','Beam S')

end
