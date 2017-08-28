%% ����RUE ����ռ��С���û���Ӱ��
clc,clear;

% ����20��
point_num = 8;
randomTotalTimeEnd720Hz = zeros(1,point_num);  % ���������ʱ��;
RateBaseTotalTimeEnd720Hz = zeros(1,point_num); % ������Ϊ��׼����ʱ��
RateBufferBaseTotalTimeEnd720Hz = zeros(1,point_num); % �����ʡ�����ռ�Ϊ��׼����ʱ��

randomTotalLossEnd720Hz = zeros(1,point_num); % �������ʧ��ʱ��
RateBaseTotalLossEnd720Hz = zeros(1,point_num); % ������Ϊ��׼ʧ��ʱ��
RateBufferBaseTotalLossEnd720Hz = zeros(1,point_num); % �����ʡ�����ռ�Ϊ��׼��ʧ��ʱ��


randomTotalTimeEnd100Hz = zeros(1,point_num);  % ���������ʱ��;
RateBaseTotalTimeEnd100Hz = zeros(1,point_num); % ������Ϊ��׼����ʱ��
RateBufferBaseTotalTimeEnd100Hz = zeros(1,point_num); % �����ʡ�����ռ�Ϊ��׼����ʱ��

randomTotalLossEnd100Hz = zeros(1,point_num); % �������ʧ��ʱ��
RateBaseTotalLossEnd100Hz = zeros(1,point_num); % ������Ϊ��׼ʧ��ʱ��
RateBufferBaseTotalLossEnd100Hz = zeros(1,point_num); % �����ʡ�����ռ�Ϊ��׼��ʧ��ʱ��

for time=1:1:5
    %% ����ÿ�η����õĻ������� 
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
    Path_loss_exponent = 4; % ·��ϵ�� ��
    K = 0.01; % ·���� K
    Power_UE = 24; % �û��ķ��书�� 24dbm    10^(24/10) mw;
    Energy_loss_factor = 0.3; %��ص�������ʧ30%
    T_slot = 1; % һ��ʱ϶ʱ�� 1us   1��=1000����(ms)1��=1000000 ΢��(��s)
    SINR_require = 20; % ��Ҫ���������� 20db

    Num_CUE = 20; % CUE�û�20��
    Num_DUE = 10; % DUE�û�10��

    OperaVol = 4; % �м̽ڵ㹤����ѹ4V
    % �ɱ����
    bandwidth = 720000; % ����720 000 hz = 720khz
    R_min = bandwidth*log2(1 + 10);

    % �����м̽ڵ������ ������Ϊ�������� 0-2000 mAh ����ֲ�
    % �������нڵ㹤��ʱ��
    RelayEnery = 2000*rand(1,10);
    t_residue = zeros(1,10);
    for i = 1:1:10
        t_residue(i) = 1000000 * judgeRelayWorkTime(RelayEnery(i),Power_UE,Energy_loss_factor,OperaVol); % us��
    end
    
    local = 1;
    %% ÿ�η���仯����
    for Relay_buffer = 1:1:8; % ����Ĵ�С1��8
        % �ɱ����
        bandwidth = 720000; % ����720 000 hz = 720khz
        %����һ ���ѡ����· 
        [tranTime,data_save_time,m1_random,m2_random] =      randomSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        %����������
        [RateBaseTime,RateBufferBaseTime,m1_RateBase,m2_RateBase,m1_RateBufferBase,m2_RateBufferBase] = comTwoWay(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        
        randomTotalTimeEnd720Hz(local) = randomTotalTimeEnd720Hz(local) + tranTime;  % ���������ʱ��;
        RateBaseTotalTimeEnd720Hz(local) = RateBaseTotalTimeEnd720Hz(local) + RateBaseTime ; % ������Ϊ��׼����ʱ��
        RateBufferBaseTotalTimeEnd720Hz(local) = RateBufferBaseTotalTimeEnd720Hz(local) + RateBufferBaseTime; % �����ʡ�����ռ�Ϊ��׼����ʱ��

        randomTotalLossEnd720Hz(local) =  randomTotalLossEnd720Hz(local) + m1_random + m2_random; % �������ʧ��ʱ��
        RateBaseTotalLossEnd720Hz(local) = RateBaseTotalLossEnd720Hz(local) + m1_RateBase + m2_RateBase; % ������Ϊ��׼ʧ��ʱ��
        RateBufferBaseTotalLossEnd720Hz(local) = RateBufferBaseTotalLossEnd720Hz(local) + m1_RateBufferBase + m2_RateBufferBase; % �����ʡ�����ռ�Ϊ��׼��ʧ��ʱ��
       
        % �ɱ����
        bandwidth = 540000; % ����540 000 hz = 540khz
        %����һ ���ѡ����· 
        [tranTime,data_save_time,m1_random,m2_random] = randomSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        %����������
        [RateBaseTime,RateBufferBaseTime,m1_RateBase,m2_RateBase,m1_RateBufferBase,m2_RateBufferBase] = comTwoWay(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        
        randomTotalTimeEnd100Hz(local) = randomTotalTimeEnd100Hz(local) + tranTime;  % ���������ʱ��;
        RateBaseTotalTimeEnd100Hz(local) = RateBaseTotalTimeEnd100Hz(local) + RateBaseTime ; % ������Ϊ��׼����ʱ��
        RateBufferBaseTotalTimeEnd100Hz(local) = RateBufferBaseTotalTimeEnd100Hz(local) + RateBufferBaseTime; % �����ʡ�����ռ�Ϊ��׼����ʱ��

        randomTotalLossEnd100Hz(local) =  randomTotalLossEnd100Hz(local) + m1_random + m2_random; % �������ʧ��ʱ��
        RateBaseTotalLossEnd100Hz(local) = RateBaseTotalLossEnd100Hz(local) + m1_RateBase + m2_RateBase; % ������Ϊ��׼ʧ��ʱ��
        RateBufferBaseTotalLossEnd100Hz(local) = RateBufferBaseTotalLossEnd100Hz(local) + m1_RateBufferBase + m2_RateBufferBase; % �����ʡ�����ռ�Ϊ��׼��ʧ��ʱ��
       
                
        local = local + 1
    end
    
    time
end

ww = 1:1:8;

%����ʱ��ͼ
randomTotalTimeEnd720Hz = randomTotalTimeEnd720Hz/time;
RateBaseTotalTimeEnd720Hz = RateBaseTotalTimeEnd720Hz/time;
RateBufferBaseTotalTimeEnd720Hz = RateBufferBaseTotalTimeEnd720Hz/time;

randomTotalTimeEnd100Hz = randomTotalTimeEnd100Hz/time;
RateBaseTotalTimeEnd100Hz = RateBaseTotalTimeEnd100Hz/time;
RateBufferBaseTotalTimeEnd100Hz = RateBufferBaseTotalTimeEnd100Hz/time;

figure(2);
plot(ww,log(randomTotalTimeEnd720Hz),'c^--'); hold on;
plot(ww,log(RateBaseTotalTimeEnd720Hz),'mo--'); hold on;
plot(ww,log(RateBufferBaseTotalTimeEnd720Hz),'rs--'); hold on;

plot(ww,log(randomTotalTimeEnd100Hz),'c^-'); hold on;
plot(ww,log(RateBaseTotalTimeEnd100Hz),'mo-'); hold on;
plot(ww,log(RateBufferBaseTotalTimeEnd100Hz),'rs-'); hold on;

% �ɹ���ͼ
randomTotalLossEnd720Hz = randomTotalLossEnd720Hz/time;
RateBaseTotalLossEnd720Hz = RateBaseTotalLossEnd720Hz/time;
RateBufferBaseTotalLossEnd720Hz = RateBufferBaseTotalLossEnd720Hz/time;


randomSucc = (randomTotalTimeEnd720Hz - randomTotalLossEnd720Hz)./randomTotalTimeEnd720Hz;
RateBaseSucc = (RateBaseTotalTimeEnd720Hz - RateBaseTotalLossEnd720Hz)./RateBaseTotalTimeEnd720Hz;
RateBufferBaseSucc = (RateBufferBaseTotalTimeEnd720Hz - RateBufferBaseTotalLossEnd720Hz)./RateBufferBaseTotalTimeEnd720Hz;

figure(3);
plot(ww,randomSucc,'c^--'); hold on;
plot(ww,RateBaseSucc,'mo--'); hold on;
plot(ww,RateBufferBaseSucc,'rs--'); hold on;

