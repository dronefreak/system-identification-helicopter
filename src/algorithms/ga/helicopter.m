function[] = helicopter()
% Example Message Types
% 
%     TIME - Time stamp
%     ATT - Vehicle attitude
%     ATSP - Vehicle attitude setpoint
%     IMU - IMU sensors
%     SENS - Other sensors
%     LPOS - Local position estimate
%     LPSP - Local position setpoint
%     GPS - GPS position
%     ATTC - Attitude controls (actuator_0 output)
%     STAT - Vehicle state
%     RC - RC input channels
%     OUT0 - Actuator_0 output
%     AIRS - Airspeed
%     ARSP - Attitude rate setpoint
%     FLOW - Optical flow
%     GPOS - Global position estimate
%     GPSP - Global position setpoint
%     ESC - ESC state
%     GVSP - Global velocity setpoint

%creates input-output variable
%IN = 1lateral 2longitudinal 3pedal 4collective
%OUT = 1roll 2pitch 3yaw 4rollrate 5pitchrate 6yawrate 7lat_vel 8lon_vel 9coll_velocity


IN(:,1) = RC(:,2);
IN(:,2) = RC(:,3);
IN(:,3) = RC(:,7);
IN(:,4) = RC(:,5);

subplot(3,1,1);
plot(IN);

OUT(:,1) = ATT(:,6);
OUT(:,2) = ATT(:,7);
OUT(:,3) = ATT(:,8);
OUT(:,4) = ATT(:,9);
OUT(:,5) = ATT(:,10);
OUT(:,6) = ATT(:,11);

subplot(3,1,2);
plot(OUT)

OUTVEL(:,1) = LPOS(:,7);
OUTVEL(:,2) = LPOS(:,8);
OUTVEL(:,3) = LPOS(:,9);

subplot(3,1,3);
plot(OUTVEL);
end