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

%load('?E:\projects\helicopter\flight data\497185.mat')

timeus = (TIME(:,2));
timeus = timeus - (min(TIME(:,2)));

IN(:,1) = RC(:,2);
IN(:,2) = RC(:,3);
IN(:,3) = RC(:,7);
IN(:,4) = RC(:,5);

subplot(4,2,1);
plot(RC(:,1),IN);
title('Lateral, Longitudinal, Pedal, Collective','FontSize',8)

OUT(:,1) = ATT(:,6);
OUT(:,2) = ATT(:,7);
OUT(:,3) = ATT(:,8);
OUT(:,4) = ATT(:,9);
OUT(:,5) = ATT(:,10);
OUT(:,6) = ATT(:,11);

subplot(3,2,3);
plot(ATT(:,1),OUT);
title('Roll, Pitch, Heading, Rollrate, Pitchrate, Yawrate','FontSize',8)

OUTVEL(:,1) = LPOS(:,7);
OUTVEL(:,2) = LPOS(:,8);
OUTVEL(:,3) = LPOS(:,9);

subplot(3,2,5);
plot(LPOS(:,1),OUTVEL);
title('Vx Vy Vz','FontSize',8)

subplot(3,2,[2,4,6])
plot3(LPOS(:,2),LPOS(:,3),LPOS(:,4));
title('Position Locus','FontSize',8)

comet3(LPOS(:,2),LPOS(:,3),LPOS(:,4))

% subplot(2,1,1);
% plot(LPOS(:,1),OUTVEL);
% subplot(2,1,2);
% plot(RC(:,1),IN);


filename = max(timeus);
fname = sprintf('HeliFD%d.mat', filename);
save(fname)