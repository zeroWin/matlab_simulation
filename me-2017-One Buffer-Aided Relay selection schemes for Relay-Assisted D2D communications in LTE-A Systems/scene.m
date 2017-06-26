%% ����ͼ����
% ��վ���Ƿ�ΧԲ
alpha = 0:pi/50:2*pi;
R_BS = 500;
x_BS = 0;
y_BS = 0;

x_cir_BS = R_BS*cos(alpha) + x_BS;
y_cir_BS = R_BS*sin(alpha) + x_BS;
plot(x_BS,y_BS,'k^','MarkerFaceColor','k');
hold on;
plot(x_cir_BS,y_cir_BS,'c--');
axis equal;
hold on;

% D2D ͨ�����긲�Ƿ�Χ�����跢�䷽��0,300�������շ���0,480��������180m
x_S = 0;
y_S = 300;
x_D = 0;
y_D = 480;


R_D2D = sqrt((x_D - x_S)^2 +(y_D - y_S)^2 );
x_cir_D2D =  R_D2D*cos(alpha)  + x_S;
y_cir_D2D =  R_D2D*sin(alpha)  + y_S;
plot(x_cir_D2D,y_cir_D2D,'r--');
axis equal;
hold on;

plot(x_S,y_S,'ro','MarkerFaceColor','r');
hold on;
plot(x_D,y_D,'rs','MarkerFaceColor','r');
hold on;
% �������CUE �û� 20��
x = 1000*rand(1,10000) - 500;
y = 1000*rand(1,10000) - 500;
x_CUE_temp = x(x.^2+y.^2<500*500);
y_CUE_temp = y(x.^2+y.^2<500*500);
x_CUE = x_CUE_temp(1:20);
y_CUE = y_CUE_temp(1:20);
scatter(x_CUE,y_CUE,'bx');
hold on;

% �������DUE �û� 10��
x = 2*R_D2D*rand(1,10000) - R_D2D;
y = 2*R_D2D*rand(1,10000) - R_D2D;
x_DUE_temp = x(x.^2+y.^2<R_D2D*R_D2D);
y_DUE_temp = y(x.^2+y.^2<R_D2D*R_D2D);
x_DUE = x_DUE_temp(1:10) + x_S;
y_DUE = y_DUE_temp(1:10) + y_S;
scatter(x_DUE,y_DUE,'mh');
hold on;
axis([-500 500 -500 500]);

%% ��������
M_packet = 50;  %������ 50
M_packet_length = 1024; % ����С 1024bit
Noise_density = -174; % �����������ܶ� -174dbm/hz    dbm = 10*lg(mw)
Path_loss_exponent = 4; % ·��ϵ��
Path_loss_constant = 0.01; % ·����
Power_UE = 24; % �û��ķ��书�� 24dbm    10^(24/10) mw;
Energy_loss_factor = 0.3; %��ص�������ʧ30%
T_slot = 1; % һ��ʱ϶ʱ�� 1us


%% ����һ ���ѡ����· 
array_fastFading =  exprnd(1);  % ��ֵΪ1��ָ���ֲ�
array_slowFading = lognrnd(0,8); %��ֵΪ0������Ϊ8db����̬�����ֲ�


%% ������ ������Ϊ׼��ѡ����·


%% ������ ͬʱ�������ʣ�����������ռ�ѡ����·



