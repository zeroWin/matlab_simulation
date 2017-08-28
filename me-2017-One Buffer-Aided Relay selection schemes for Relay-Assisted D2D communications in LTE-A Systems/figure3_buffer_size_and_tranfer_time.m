%% 仿真RUE 缓存空间大小对用户的影响
clc,clear;

% 仿真20次
point_num = 8;
randomTotalTimeEnd720Hz = zeros(1,point_num);  % 随机方法总时间;
RateBaseTotalTimeEnd720Hz = zeros(1,point_num); % 以速率为基准的总时间
RateBufferBaseTotalTimeEnd720Hz = zeros(1,point_num); % 以速率、缓存空间为基准的总时间

randomTotalLossEnd720Hz = zeros(1,point_num); % 随机方法失败时间
RateBaseTotalLossEnd720Hz = zeros(1,point_num); % 以速率为基准失败时间
RateBufferBaseTotalLossEnd720Hz = zeros(1,point_num); % 以速率、缓存空间为基准的失败时间


randomTotalTimeEnd100Hz = zeros(1,point_num);  % 随机方法总时间;
RateBaseTotalTimeEnd100Hz = zeros(1,point_num); % 以速率为基准的总时间
RateBufferBaseTotalTimeEnd100Hz = zeros(1,point_num); % 以速率、缓存空间为基准的总时间

randomTotalLossEnd100Hz = zeros(1,point_num); % 随机方法失败时间
RateBaseTotalLossEnd100Hz = zeros(1,point_num); % 以速率为基准失败时间
RateBufferBaseTotalLossEnd100Hz = zeros(1,point_num); % 以速率、缓存空间为基准的失败时间

for time=1:1:5
    %% 生成每次仿真用的基本参数 
    % 基站覆盖范围圆
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

    % D2D 通信坐标覆盖范围，假设发射方在0,300处，接收方在0,480处，距离180m
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
    % 随机产生CUE 用户 20个
    x = 1000*rand(1,10000) - 500;
    y = 1000*rand(1,10000) - 500;
    x_CUE_temp = x(x.^2+y.^2<500*500);
    y_CUE_temp = y(x.^2+y.^2<500*500);
    x_CUE = x_CUE_temp(1:20);
    y_CUE = y_CUE_temp(1:20);
    scatter(x_CUE,y_CUE,'bx');
    hold on;

    % 随机产生DUE 用户 10个
    x = 2*R_D2D*rand(1,10000) - R_D2D;
    y = 2*R_D2D*rand(1,10000) - R_D2D;
    x_DUE_temp = x(x.^2+y.^2<R_D2D*R_D2D);
    y_DUE_temp = y(x.^2+y.^2<R_D2D*R_D2D);
    x_DUE = x_DUE_temp(1:10) + x_S;
    y_DUE = y_DUE_temp(1:10) + y_S;
    scatter(x_DUE,y_DUE,'mh');
    hold on;
    axis([-500 500 -500 500]);

    %% 基本参数
    M_packet = 50;  %包个数 50
    M_packet_length = 1024; % 包大小 1024bit
    Noise_density = -174; % 噪声功率谱密度 -174dbm/hz    dbm = 10*lg(mw)
    Path_loss_exponent = 4; % 路损系数 α
    K = 0.01; % 路损常数 K
    Power_UE = 24; % 用户的发射功率 24dbm    10^(24/10) mw;
    Energy_loss_factor = 0.3; %电池的能量损失30%
    T_slot = 1; % 一个时隙时常 1us   1秒=1000毫秒(ms)1秒=1000000 微秒(μs)
    SINR_require = 20; % 需要的最低信噪比 20db

    Num_CUE = 20; % CUE用户20个
    Num_DUE = 10; % DUE用户10个

    OperaVol = 4; % 中继节点工作电压4V
    % 可变参数
    bandwidth = 720000; % 带宽720 000 hz = 720khz
    R_min = bandwidth*log2(1 + 10);

    % 生成中继节点的能量 暂且认为能量无限 0-2000 mAh 随机分布
    % 计算所有节点工作时间
    RelayEnery = 2000*rand(1,10);
    t_residue = zeros(1,10);
    for i = 1:1:10
        t_residue(i) = 1000000 * judgeRelayWorkTime(RelayEnery(i),Power_UE,Energy_loss_factor,OperaVol); % us数
    end
    
    local = 1;
    %% 每次仿真变化参数
    for Relay_buffer = 1:1:8; % 缓存的大小1到8
        % 可变参数
        bandwidth = 720000; % 带宽720 000 hz = 720khz
        %方法一 随机选择链路 
        [tranTime,data_save_time,m1_random,m2_random] =      randomSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        %方法二、三
        [RateBaseTime,RateBufferBaseTime,m1_RateBase,m2_RateBase,m1_RateBufferBase,m2_RateBufferBase] = comTwoWay(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        
        randomTotalTimeEnd720Hz(local) = randomTotalTimeEnd720Hz(local) + tranTime;  % 随机方法总时间;
        RateBaseTotalTimeEnd720Hz(local) = RateBaseTotalTimeEnd720Hz(local) + RateBaseTime ; % 以速率为基准的总时间
        RateBufferBaseTotalTimeEnd720Hz(local) = RateBufferBaseTotalTimeEnd720Hz(local) + RateBufferBaseTime; % 以速率、缓存空间为基准的总时间

        randomTotalLossEnd720Hz(local) =  randomTotalLossEnd720Hz(local) + m1_random + m2_random; % 随机方法失败时间
        RateBaseTotalLossEnd720Hz(local) = RateBaseTotalLossEnd720Hz(local) + m1_RateBase + m2_RateBase; % 以速率为基准失败时间
        RateBufferBaseTotalLossEnd720Hz(local) = RateBufferBaseTotalLossEnd720Hz(local) + m1_RateBufferBase + m2_RateBufferBase; % 以速率、缓存空间为基准的失败时间
       
        % 可变参数
        bandwidth = 540000; % 带宽540 000 hz = 540khz
        %方法一 随机选择链路 
        [tranTime,data_save_time,m1_random,m2_random] = randomSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        %方法二、三
        [RateBaseTime,RateBufferBaseTime,m1_RateBase,m2_RateBase,m1_RateBufferBase,m2_RateBufferBase] = comTwoWay(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        
        randomTotalTimeEnd100Hz(local) = randomTotalTimeEnd100Hz(local) + tranTime;  % 随机方法总时间;
        RateBaseTotalTimeEnd100Hz(local) = RateBaseTotalTimeEnd100Hz(local) + RateBaseTime ; % 以速率为基准的总时间
        RateBufferBaseTotalTimeEnd100Hz(local) = RateBufferBaseTotalTimeEnd100Hz(local) + RateBufferBaseTime; % 以速率、缓存空间为基准的总时间

        randomTotalLossEnd100Hz(local) =  randomTotalLossEnd100Hz(local) + m1_random + m2_random; % 随机方法失败时间
        RateBaseTotalLossEnd100Hz(local) = RateBaseTotalLossEnd100Hz(local) + m1_RateBase + m2_RateBase; % 以速率为基准失败时间
        RateBufferBaseTotalLossEnd100Hz(local) = RateBufferBaseTotalLossEnd100Hz(local) + m1_RateBufferBase + m2_RateBufferBase; % 以速率、缓存空间为基准的失败时间
       
                
        local = local + 1
    end
    
    time
end

ww = 1:1:8;

%传输时间图
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

% 成功率图
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

