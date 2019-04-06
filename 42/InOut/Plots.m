R2D = 180/pi;
D2R = pi/180;

fileID = fopen('time.42','r');
formatSpec = '%f';
sizeA = [1 Inf];
time = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('Cmd.42','r');
formatSpec = '%f';
sizeA = [4 Inf];
cmd = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('acqbn.42','r');
formatSpec = '%f';
sizeA = [4 Inf];
acqbn = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('Cmd2.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
cmd2 = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('PosN.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
PosN = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('WhlCmd.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
WhlCmd = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('MtqCmd.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
MtqCmd = fscanf(fileID, formatSpec, sizeA);

% fileID = fopen('PosNT1.42','r');
% formatSpec = '%f';
% sizeA = [3 Inf];
% PosNT1 = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('PosNT2.42','r');
% formatSpec = '%f';
% sizeA = [3 Inf];
% PosNT2 = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('VelN.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
VelN = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('wbn.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
wbn = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('qbn.42','r');
formatSpec = '%f';
sizeA = [4 Inf];
qbn = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('wbn2.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
wbn2 = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('qbn2.42','r');
formatSpec = '%f';
sizeA = [4 Inf];
qbn2 = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('wbn3.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
wbn3 = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('qbn3.42','r');
formatSpec = '%f';
sizeA = [4 Inf];
qbn3 = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('Hvn.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
Hvn = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('RPY.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
RPY = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('Gyro.42','r');
% formatSpec = '%f';
% sizeA = [18 Inf];
% Gyro = fscanf(fileID, formatSpec, sizeA);
% Gyro = Gyro * 180 / pi;
% gyro_wn = Gyro(1:3,:);
% gyro_TrueRate = Gyro(4:6,:);
% gyro_Bias = Gyro(7:9, :);
% gyro_Angle = Gyro(10:12, :);
% gyro_MeasRate = Gyro(13:15, :);
% gyro_ACwbn = Gyro(16:18, :);

fileID = fopen('alt.42','r');
formatSpec = '%f';
sizeA = [1 Inf];
alt = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('EclipseFlag.42','r');
formatSpec = '%d';
sizeA = [1 Inf];
EclipseFlag = fscanf(fileID, formatSpec, sizeA);

fileID = fopen('PosW.42','r');
formatSpec = '%f';
sizeA = [3 Inf];
PosWR = fscanf(fileID, formatSpec, sizeA);

% fileID = fopen('PosWT1.42','r');
% formatSpec = '%f';
% sizeA = [3 Inf];
% PosWT1 = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('PosWT2.42','r');
% formatSpec = '%f';
% sizeA = [3 Inf];
% PosWT2 = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('WGS84T1.42','r');
% formatSpec = '%f';
% sizeA = [3 Inf];
% WGS84T1 = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('WGS84T2.42','r');
% formatSpec = '%f';
% sizeA = [3 Inf];
% WGS84T2 = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('WGS84R.42','r');
% formatSpec = '%f';
% sizeA = [3 Inf];
% WGS84R = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('SP1.42','r');
% formatSpec = '%f';
% sizeA = [7 Inf];
% SP1 = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('SP2.42','r');
% formatSpec = '%f';
% sizeA = [7 Inf];
% SP2 = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('Swath.42','r');
% formatSpec = '%f';
% sizeA = [2 Inf];
% SwathArea = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('PathLength1.42','r');
% formatSpec = '%f';
% sizeA = [1 Inf];
% PathLength1 = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('PathLength2.42','r');
% formatSpec = '%f';
% sizeA = [1 Inf];
% PathLength2 = fscanf(fileID, formatSpec, sizeA);
% 
% fileID = fopen('ValidSp.42','r');
% formatSpec = '%f';
% sizeA = [5 Inf];
% ValidSp = fscanf(fileID, formatSpec, sizeA);

% figure
% yyaxis left
% plot(time, SwathArea(1,:)/1000000);
% title('Swath Area & Width');
% xlabel('time[sec]');
% ylabel('Area[km^2]');
% hold on
% yyaxis right
% plot(time, SwathArea(2,:)/1000);
% ylabel('Width[km]');
% hold off

% figure
% plot(time,alt/1000)
% title('Altitude of SC')
% xlabel('time[sec]')
% ylabel('altitude[km]')
% grid on

% [row, col] = find(~PosNT2);
% for i=1:length(row)
%     for j=1:length(col)
%         PosNT2(row(i),col(j)) = NaN;
%         WGS84T2(row(i),col(j)) = NaN;
%         SP2(row(i),col(j)) = NaN;
%         SP2(row(i)+3,col(j)) = NaN;
%     end
% end

% figure
% subplot(3,1,1);
% plot(time, PosN(1,:), time, PosNT1(1,:),'g', time, SP1(1,:), 'r');
% title('Position of Rx & MUOS & Specular Point wrt ECI Frame');
% xlabel('time[sec]');
% ylabel('x[m]');
% subplot(3,1,2);
% plot(time, PosN(2,:), time, PosNT1(2,:),'g', time, SP1(2,:), 'r');
% xlabel('time[sec]');
% ylabel('y[m]');
% subplot(3,1,3);
% plot(time, PosN(3,:), time, PosNT1(3,:),'g', time, SP1(3,:), 'r');
% xlabel('time[sec]');
% ylabel('z[m]');
% legend('Rx','Tx','SP');
% 
% figure
% subplot(3,1,1);
% plot(time, PosN(1,:), time, PosNT2(1,:),'g', time, SP2(1,:), 'r');
% title('Position of Rx & Orbcomms & Specular Point wrt ECI Frame');
% xlabel('time[sec]');
% ylabel('x[m]');
% subplot(3,1,2);
% plot(time, PosN(2,:), time, PosNT2(2,:),'g', time, SP2(2,:), 'r');
% xlabel('time[sec]');
% ylabel('y[m]');
% subplot(3,1,3);
% plot(time, PosN(3,:), time, PosNT2(3,:),'g', time, SP2(3,:), 'r');
% xlabel('time[sec]');
% ylabel('z[m]');
% legend('Rx','Tx','SP');
% 
% figure
% subplot(3,1,1);
% plot(time, WGS84R(1,:)*R2D, time, WGS84T1(1,:)*R2D,'g', time, SP1(4,:)*R2D, 'r');
% title('LLA of Rx & MUOS & Specular Point');
% xlabel('time[sec]');
% ylabel('lat[deg]');
% ylim([-55 55])
% subplot(3,1,2);
% plot(time, WGS84R(2,:)*R2D, time, WGS84T1(2,:)*R2D,'g', time, SP1(5,:)*R2D, 'r');
% xlabel('time[sec]');
% ylabel('lon[deg]');
% subplot(3,1,3);
% plot(time, WGS84R(3,:)/1000, time, WGS84T1(3,:)/1000,'g', time, SP1(6,:)/1000, 'r');
% xlabel('time[sec]');
% ylabel('alt[km]');
% legend('Rx','Tx','SP');
% 
% figure
% subplot(3,1,1);
% plot(time, WGS84R(1,:)*R2D, time, WGS84T2(1,:)*R2D,'g', time, SP2(4,:)*R2D, 'r');
% title('LLA of Rx & Orbcomms & Specular Point');
% xlabel('time[sec]');
% ylabel('lat[deg]');
% ylim([-55 55])
% subplot(3,1,2);
% plot(time, WGS84R(2,:)*R2D, time, WGS84T2(2,:)*R2D,'g', time, SP2(5,:)*R2D, 'r');
% xlabel('time[sec]');
% ylabel('lon[deg]');
% subplot(3,1,3);
% plot(time, WGS84R(3,:)/1000, time, WGS84T2(3,:)/1000,'g', time, SP2(6,:)/1000, 'r');
% xlabel('time[sec]');
% ylabel('alt[km]');
% legend('Rx','Tx','SP');

% for i=1:length(time)
%     RE(i)=6.378145E6;
%     R(i)=sqrt(PosN(1,i)^2+PosN(2,i)^2+PosN(3,i)^2)-RE(i);
% end
% 
% plot(time,R)



% figure
% subplot(3,1,1);
% plot(time,PosN(1,:));
% title('Position of SC wrt ECI Frame');
% xlabel('time[sec]');
% ylabel('x[m]');
% subplot(3,1,2);
% plot(time,PosN(2,:));
% xlabel('time[sec]');
% ylabel('y[m]');
% subplot(3,1,3);
% plot(time,PosN(3,:));
% xlabel('time[sec]');
% ylabel('z[m]');
% 
% figure
% subplot(3,1,1);
% plot(time,VelN(1,:));
% title('Velocity of SC wrt ECI Frame');
% xlabel('time[sec]');
% ylabel('v_x[m/s]');
% subplot(3,1,2);
% plot(time,VelN(2,:));
% xlabel('time[sec]');
% ylabel('v_y[m/s]');
% subplot(3,1,3);
% plot(time,VelN(3,:));
% xlabel('time[sec]');
% ylabel('v_z[m/s]');

figure
subplot(3,1,1);
plot(time,Hvn(1,:));
title('Total Angular Momentum of SC wrt Inertial Frame');
xlabel('time[sec]');
ylabel('h_x[Nms]');
subplot(3,1,2);
plot(time,Hvn(2,:));
xlabel('time[sec]');
ylabel('h_y[Nms]');
subplot(3,1,3);
plot(time,Hvn(3,:));
xlabel('time[sec]');
ylabel('h_z[Nms]');

figure
subplot(3,1,1)
plot(time,wbn(1,:)*R2D , time , cmd2(1,:)*R2D, '--')
title('Angular Velocity of Body 1 wrt Body 1 frame')
xlabel('time[sec]')
ylabel('w_x[deg/s]')
subplot(3,1,2)
plot(time,wbn(2,:)*R2D , time , cmd2(2,:)*R2D, '--')
xlabel('time[sec]')
ylabel('w_y[deg/s]')
subplot(3,1,3)
plot(time,wbn(3,:)*R2D , time , cmd2(3,:)*R2D, '--')
xlabel('time[sec]')
ylabel('w_z[deg/s]')
legend('True', 'Command')

figure
subplot(3,1,1);
plot(time,wbn(1,:)-cmd2(1,:));
title('Angular Velocity Error of Body 1 wrt Body 1 frame');
xlabel('time[sec]');
ylabel('w_x[rad/s]');
subplot(3,1,2);
plot(time,wbn(2,:)-cmd2(2,:));
xlabel('time[sec]');
ylabel('w_y[rad/s]');
subplot(3,1,3);
plot(time,wbn(3,:)-cmd2(3,:));
xlabel('time[sec]');
ylabel('w_z[rad/s]');

% figure
% subplot(3,1,1);
% plot(time,wbn(1,:) , time , wbn2(1,:), '--', time , wbn3(1,:), '-.');
% title('Angular Velocity of SC wrt Body 1, 2, 3');
% xlabel('time[sec]');
% ylabel('w_x[rad/s]');
% subplot(3,1,2);
% plot(time,wbn(2,:) , time , wbn2(2,:), '--', time , wbn3(2,:), '-.');
% xlabel('time[sec]');
% ylabel('w_y[rad/s]');
% subplot(3,1,3);
% plot(time,wbn(3,:) , time , wbn2(3,:), '--', time , wbn3(3,:), '-.');
% xlabel('time[sec]');
% ylabel('w_z[rad/s]');
% legend('Body 1', 'Body 2', 'Body 3');
% 
% figure
% subplot(4,1,1);
% plot(time,qbn(1,:), time, qbn2(1,:), '--', time, qbn3(1,:), '-.');
% title('Quaternion of SC wrt Body 1, 2, 3');
% xlabel('time[sec]');
% ylabel('q_1');
% subplot(4,1,2);
% plot(time,qbn(2,:), time, qbn2(2,:), '--', time, qbn3(2,:), '-.');
% xlabel('time[sec]');
% ylabel('q_2');
% subplot(4,1,3);
% plot(time,qbn(3,:), time, qbn2(3,:), '--', time, qbn3(3,:), '-.');
% xlabel('time[sec]');
% ylabel('q_3');
% subplot(4,1,4);
% plot(time,qbn(4,:), time, qbn2(4,:), '--', time, qbn3(4,:), '-.');
% xlabel('time[sec]');
% ylabel('q_4');
% legend('Body 1', 'Body 2', 'Body 3');

% figure
% subplot(3,1,1);
% plot(time,wbn(1,:)-wbn2(1,:), time,wbn(1,:)-wbn3(1,:), '--');
% title('Angular Velocity Error');
% xlabel('time[sec]');
% ylabel('w_x[rad/s]');
% subplot(3,1,2);
% plot(time,wbn(2,:)-wbn2(2,:), time,wbn(2,:)-wbn3(2,:), '--');
% xlabel('time[sec]');
% ylabel('w_y[rad/s]');
% subplot(3,1,3);
% plot(time,wbn(3,:)-wbn2(3,:), time,wbn(3,:)-wbn3(3,:), '--');
% xlabel('time[sec]');
% ylabel('w_z[rad/s]');
% legend('Body 1 - Body 2', 'Body 1 - Body 3');
% 
% figure
% subplot(4,1,1);
% plot(time,qbn(1,:)-qbn2(1,:));
% title('Quaternion Difference btw Body 1 & 2');
% xlabel('time[sec]');
% ylabel('q_1');
% subplot(4,1,2);
% plot(time,qbn(2,:)-qbn2(2,:));
% xlabel('time[sec]');
% ylabel('q_2');
% subplot(4,1,3);
% plot(time,qbn(3,:)-qbn2(3,:));
% xlabel('time[sec]');
% ylabel('q_3');
% subplot(4,1,4);
% plot(time,qbn(4,:)-qbn2(4,:));
% xlabel('time[sec]');
% ylabel('q_4');

figure
subplot(3,1,1);
plot(time,RPY(1,:));
title('Euler Angles of Body 1 wrt LVLH Frame');
xlabel('time[sec]');
ylabel('roll[deg]');
subplot(3,1,2);
plot(time,RPY(2,:));
xlabel('time[sec]');
ylabel('pitch[deg]');
subplot(3,1,3);
plot(time,RPY(3,:));
xlabel('time[sec]');
ylabel('yaw[deg]');

figure
subplot(3,1,1);
plot(time,MtqCmd(1,:));
title('MTQ torque command');
xlabel('time[sec]');
ylabel('M');
subplot(3,1,2);
plot(time,MtqCmd(2,:));
xlabel('time[sec]');
ylabel('M');
subplot(3,1,3);
plot(time,MtqCmd(3,:));
xlabel('time[sec]');
ylabel('M');

figure
subplot(3,1,1);
plot(time,WhlCmd(1,:));
title('Reaction wheel torque command');
xlabel('time[sec]');
ylabel('M');
subplot(3,1,2);
plot(time,WhlCmd(2,:));
xlabel('time[sec]');
ylabel('M');
subplot(3,1,3);
plot(time,WhlCmd(3,:));
xlabel('time[sec]');
ylabel('M');

figure
subplot(4,1,1);
plot(time,cmd(1,:), time,qbn(1,:), time,acqbn(1,:));
title('Quaternion command');
xlabel('time[sec]');
ylabel('M');
subplot(4,1,2);
plot(time,cmd(2,:), time,qbn(2,:), time,acqbn(2,:));
xlabel('time[sec]');
ylabel('M');
subplot(4,1,3);
plot(time,cmd(3,:), time,qbn(3,:), time,acqbn(3,:));
xlabel('time[sec]');
ylabel('M');
subplot(4,1,4);
plot(time,cmd(4,:), time,qbn(4,:), time,acqbn(4,:));
xlabel('time[sec]');
ylabel('M');
legend('cmd_qrn','qbn', 'acqbn', 'stqbn', 'stval', 'BoS');

figure
subplot(3,1,1);
plot(time,cmd2(1,:), time,wbn(1,:));
title('Reaction wheel torque command');
xlabel('time[sec]');
ylabel('M');
subplot(3,1,2);
plot(time,cmd2(2,:), time,wbn(2,:));
xlabel('time[sec]');
ylabel('M');
subplot(3,1,3);
plot(time,cmd2(3,:), time,wbn(3,:));
xlabel('time[sec]');
ylabel('M');


% figure
% plot(time,PathLength1/1000)
% title('Signal Path Length from MUOS')
% xlabel('time[sec]')
% ylabel('distance[km]')
% grid on
% 
% figure
% plot(time,PathLength2/1000)
% title('Signal Path Distance from Orbcomm')
% xlabel('time[sec]')
% ylabel('distance[km]')
% grid on
% figure
% subplot(3,1,1);
% plot(time,gyro_wn(1,:), time,gyro_ACwbn(1,:), time, gyro_wn(1,:) - gyro_ACwbn(1,:));
% title('Gyro Angluar Velocity');
% xlabel('time[sec]');
% ylabel('x[deg/sec]');
% subplot(3,1,2);
% plot(time,gyro_wn(2,:), time,gyro_ACwbn(2,:), time, gyro_wn(2,:) - gyro_ACwbn(2,:));
% xlabel('time[sec]');
% ylabel('y[deg/sec]');
% subplot(3,1,3);
% plot(time,gyro_wn(3,:), time,gyro_ACwbn(3,:), time, gyro_wn(3,:) - gyro_ACwbn(3,:));
% xlabel('time[sec]');
% ylabel('z[deg/sec]');

% figure
% subplot(3,1,1);
% plot(time,gyro_Bias(1,:), time,gyro_MeasRate(1,:), time, gyro_TrueRate(1,:));
% title('Gyro Rate');
% xlabel('time[sec]');
% ylabel('x[deg/sec]');
% subplot(3,1,2);
% plot(time,gyro_Bias(2,:), time,gyro_MeasRate(2,:), time, gyro_TrueRate(2,:));
% xlabel('time[sec]');
% ylabel('y[deg/sec]');
% subplot(3,1,3);
% plot(time,gyro_Bias(3,:), time,gyro_MeasRate(3,:), time, gyro_TrueRate(3,:));
% xlabel('time[sec]');
% ylabel('z[deg/sec]');
% 
% figure
% subplot(3,1,1);
% plot(time,gyro_Angle(1,:));
% title('Gyro Angle');
% xlabel('time[sec]');
% ylabel('x[deg]');
% subplot(3,1,2);
% plot(time,gyro_Angle(2,:));
% xlabel('time[sec]');
% ylabel('y[deg]');
% subplot(3,1,3);
% plot(time,gyro_Angle(3,:));
% xlabel('time[sec]');
% ylabel('z[deg]');
% 
% figure
% subplot(3,1,1);
% plot(time,gyro_ACwbn(1,:)-gyro_MeasRate(1,:), time, gyro_wn(1,:) - gyro_MeasRate(1,:));
% title('Gyro Angluar Velocity');
% xlabel('time[sec]');
% ylabel('x[deg/sec]');
% subplot(3,1,2);
% plot(time,gyro_ACwbn(2,:)-gyro_MeasRate(2,:), time, gyro_wn(2,:) - gyro_MeasRate(2,:));
% xlabel('time[sec]');
% ylabel('y[deg/sec]');
% subplot(3,1,3);
% plot(time,gyro_ACwbn(3,:)-gyro_MeasRate(3,:), time, gyro_wn(3,:) - gyro_MeasRate(3,:));
% xlabel('time[sec]');
% ylabel('z[deg/sec]');