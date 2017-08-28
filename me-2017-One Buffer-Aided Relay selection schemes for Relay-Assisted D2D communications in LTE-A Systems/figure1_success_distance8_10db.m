%% ����S->D֮��ľ����ͨ��ʱ���Ӱ��
clc,clear;
%% ����ͼ����
% ��վ���Ƿ�ΧԲ
Num_CUE = 20; % CUE�û�20��
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

% �������CUE �û� 20��
x = 1000*rand(1,10000) - 500;
y = 1000*rand(1,10000) - 500;
x_CUE_temp = x(x.^2+y.^2<500*500);
y_CUE_temp = y(x.^2+y.^2<500*500);
x_CUE = x_CUE_temp(1:Num_CUE);
y_CUE = y_CUE_temp(1:Num_CUE);
scatter(x_CUE,y_CUE,'bx');
hold on;

%% ��������
M_packet = 50;  %������ 50
M_packet_length = 1024; % ����С 1024bit
Relay_buffer = 5; % �ܴ洢������5
Noise_density = -174; % �����������ܶ� -174dbm/hz    dbm = 10*lg(mw)
Path_loss_exponent = 4; % ·��ϵ�� ��
K = 0.01; % ·���� K
Power_UE = 24; % �û��ķ��书�� 24dbm    10^(24/10) mw;
Energy_loss_factor = 0.3; %��ص�������ʧ30%
T_slot = 1; % һ��ʱ϶ʱ�� 1us   1��=1000����(ms)1��=1000000 ΢��(��s)
SINR_require = 8; % ��Ҫ���������� 10db


Num_DUE = 10; % DUE�û�10��
 
OperaVol = 4; % �м̽ڵ㹤����ѹ4V
% �ɱ����
bandwidth = 720000; % ����720 000 hz = 720khz
%bandwidth = 1; % ����720 000 hz = 720khz

 %% S��D֮������ͨ�ųɹ��ʵ�Ӱ��
 % D2D ͨ�����긲�Ƿ�Χ�����跢�䷽��0,250�������շ���0,470��������220m
x_S = 0;
y_S = 250;
x_D = 0;

point_num = 6;
% 8db
randomTotalTimeEnd8db = zeros(1,point_num);  % ���������ʱ��;
RateBaseTotalTimeEnd8db = zeros(1,point_num); % ������Ϊ��׼����ʱ��
RateBufferBaseTotalTimeEnd8db = zeros(1,point_num); % �����ʡ�����ռ�Ϊ��׼����ʱ��

randomTotalLossEnd8db = zeros(1,point_num); % �������ʧ��ʱ��
RateBaseTotalLossEnd8db = zeros(1,point_num); % ������Ϊ��׼ʧ��ʱ��
RateBufferBaseTotalLossEnd8db = zeros(1,point_num); % �����ʡ�����ռ�Ϊ��׼��ʧ��ʱ��

% 10db
randomTotalTimeEnd10db = zeros(1,point_num);  % ���������ʱ��;
RateBaseTotalTimeEnd10db = zeros(1,point_num); % ������Ϊ��׼����ʱ��
RateBufferBaseTotalTimeEnd10db = zeros(1,point_num); % �����ʡ�����ռ�Ϊ��׼����ʱ��

randomTotalLossEnd10db = zeros(1,point_num); % �������ʧ��ʱ��
RateBaseTotalLossEnd10db = zeros(1,point_num); % ������Ϊ��׼ʧ��ʱ��
RateBufferBaseTotalLossEnd10db = zeros(1,point_num); % �����ʡ�����ռ�Ϊ��׼��ʧ��ʱ��


for time =  1:1:10
    
    local = 1;
    for y_D = 270:40:470 % ����20-220  270-40-470
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
        % �������DUE �û� 10��
        x = 2*R_D2D*rand(1,10000) - R_D2D;
        y = 2*R_D2D*rand(1,10000) - R_D2D;
        x_DUE_temp = x(x.^2+y.^2<R_D2D*R_D2D);
        y_DUE_temp = y(x.^2+y.^2<R_D2D*R_D2D);
        x_DUE = x_DUE_temp(1:Num_DUE) + x_S;
        y_DUE = y_DUE_temp(1:Num_DUE) + y_S;
        scatter(x_DUE,y_DUE,'mh');
        hold on;
        axis([-500 500 -500 500]);


        % �����м̽ڵ������ ������Ϊ�������� 0-2000 mAh ����ֲ�
        % �������нڵ㹤��ʱ��
        RelayEnery = 2000*rand(1,Num_DUE);
        t_residue = zeros(1,Num_DUE);
        for i = 1:1:Num_DUE
            t_residue(i) = 1000000 * judgeRelayWorkTime(RelayEnery(i),Power_UE,Energy_loss_factor,OperaVol); % us��
        end
        
        % 8db
        SINR_require = 70; % ��Ҫ���������� 8db
        %����һ ���ѡ����· 
        [tranTime,data_save_time,m1_random,m2_random,t_last_random] = randomSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        %����������
        [RateBaseTime,RateBufferBaseTime,m1_RateBase,m2_RateBase,m1_RateBufferBase,m2_RateBufferBase,t_last_RataBase,t_last_RateBufferBase] = comTwoWay(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);

        randomTotalTimeEnd8db(local) =  randomTotalTimeEnd8db(local) + tranTime;
        randomTotalLossEnd8db(local) = randomTotalLossEnd8db(local)  + m1_random + m2_random;

        RateBaseTotalTimeEnd8db(local) = RateBaseTotalTimeEnd8db(local) +  RateBaseTime;
        RateBaseTotalLossEnd8db(local) =  RateBaseTotalLossEnd8db(local) +  m1_RateBase + m2_RateBase;

        RateBufferBaseTotalTimeEnd8db(local) =  RateBufferBaseTotalTimeEnd8db(local) + RateBufferBaseTime;
        RateBufferBaseTotalLossEnd8db(local) = RateBufferBaseTotalLossEnd8db(local) + m1_RateBufferBase + m2_RateBufferBase;
   
         % 10db
        SINR_require = 20; % ��Ҫ���������� 10db
        %����һ ���ѡ����· 
        [tranTime,data_save_time,m1_random,m2_random,t_last_random] = randomSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        %����������
        [RateBaseTime,RateBufferBaseTime,m1_RateBase,m2_RateBase,m1_RateBufferBase,m2_RateBufferBase,t_last_RataBase,t_last_RateBufferBase] = comTwoWay(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);

        randomTotalTimeEnd10db(local) =  randomTotalTimeEnd10db(local) + tranTime;
        randomTotalLossEnd10db(local) = randomTotalLossEnd10db(local)  + m1_random + m2_random;

        RateBaseTotalTimeEnd10db(local) = RateBaseTotalTimeEnd10db(local) +  RateBaseTime;
        RateBaseTotalLossEnd10db(local) =  RateBaseTotalLossEnd10db(local) +  m1_RateBase + m2_RateBase;

        RateBufferBaseTotalTimeEnd10db(local) =  RateBufferBaseTotalTimeEnd10db(local) + RateBufferBaseTime;
        RateBufferBaseTotalLossEnd10db(local) = RateBufferBaseTotalLossEnd10db(local) + m1_RateBufferBase + m2_RateBufferBase;
        
        local = local + 1
    end
    time
end
% % ����ͨ�ųɹ���
% randomSucc = (randomTotalTime - randomLoss)./randomTotalTime;
% RateBaseSucc = (RateBaseTotalTime - RateBaseLoss)./RateBaseTotalTime;
% RateBufferBaseSucc = (RateBufferBaseTotalTime - RateBufferBaseTotalLoss)./RateBufferBaseTotalTime;
% 8db
randomTotalTimeEnd8db = randomTotalTimeEnd8db/time;
RateBaseTotalTimeEnd8db = RateBaseTotalTimeEnd8db/time;
RateBufferBaseTotalTimeEnd8db = RateBufferBaseTotalTimeEnd8db/time;

randomTotalLossEnd8db = randomTotalLossEnd8db/time;
RateBaseTotalLossEnd8db = RateBaseTotalLossEnd8db/time;
RateBufferBaseTotalLossEnd8db = RateBufferBaseTotalLossEnd8db/time;

% 10db
randomTotalTimeEnd10db = randomTotalTimeEnd10db/time;
RateBaseTotalTimeEnd10db = RateBaseTotalTimeEnd10db/time;
RateBufferBaseTotalTimeEnd10db = RateBufferBaseTotalTimeEnd10db/time;

randomTotalLossEnd10db = randomTotalLossEnd10db/time;
RateBaseTotalLossEnd10db = RateBaseTotalLossEnd10db/time;
RateBufferBaseTotalLossEnd10db = RateBufferBaseTotalLossEnd10db/time;

% ��ͼ
ww = 20:40:220;
%ww = 20:20:80;
figure(2);
plot(ww,log(randomTotalTimeEnd8db),'c^--'); hold on;
plot(ww,log(RateBaseTotalTimeEnd8db),'mo--'); hold on;
plot(ww,log(RateBufferBaseTotalTimeEnd8db),'rs--'); hold on;

plot(ww,log(randomTotalTimeEnd10db),'b^-'); hold on;
plot(ww,log(RateBaseTotalTimeEnd10db),'ko-'); hold on;
plot(ww,log(RateBufferBaseTotalTimeEnd10db),'gs-'); hold on;

%8db
randomSucc8db = (randomTotalTimeEnd8db - randomTotalLossEnd8db)./randomTotalTimeEnd8db;
RateBaseSucc8db = (RateBaseTotalTimeEnd8db - RateBaseTotalLossEnd8db)./RateBaseTotalTimeEnd8db;
RateBufferBaseSucc8db = (RateBufferBaseTotalTimeEnd8db - RateBufferBaseTotalLossEnd8db)./RateBufferBaseTotalTimeEnd8db;

randomSucc10db = (randomTotalTimeEnd10db - randomTotalLossEnd10db)./randomTotalTimeEnd10db;
RateBaseSucc10db = (RateBaseTotalTimeEnd10db - RateBaseTotalLossEnd10db)./RateBaseTotalTimeEnd10db;
RateBufferBaseSucc10db = (RateBufferBaseTotalTimeEnd10db - RateBufferBaseTotalLossEnd10db)./RateBufferBaseTotalTimeEnd10db;


figure(3);
plot(ww,randomSucc8db,'c^--'); hold on;
plot(ww,RateBaseSucc8db,'mo--'); hold on;
plot(ww,RateBufferBaseSucc8db,'rs--'); hold on;

plot(ww,randomSucc10db,'b^-'); hold on;
plot(ww,RateBaseSucc10db,'ko-'); hold on;
plot(ww,RateBufferBaseSucc10db,'gs-'); hold on;
 