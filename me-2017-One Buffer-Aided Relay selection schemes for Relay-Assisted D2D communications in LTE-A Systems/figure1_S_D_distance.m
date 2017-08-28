%% 仿真S->D之间的距离对通信时间的影响
clc,clear;
%% 场景图仿真
% 基站覆盖范围圆
Num_CUE = 20; % CUE用户20个
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

% 随机产生CUE 用户 20个
x = 1000*rand(1,10000) - 500;
y = 1000*rand(1,10000) - 500;
x_CUE_temp = x(x.^2+y.^2<500*500);
y_CUE_temp = y(x.^2+y.^2<500*500);
x_CUE = x_CUE_temp(1:Num_CUE);
y_CUE = y_CUE_temp(1:Num_CUE);
scatter(x_CUE,y_CUE,'bx');
hold on;

%% 基本参数
M_packet = 50;  %包个数 50
M_packet_length = 1024; % 包大小 1024bit
Relay_buffer = 5; % 能存储包个数5
Noise_density = -174; % 噪声功率谱密度 -174dbm/hz    dbm = 10*lg(mw)
Path_loss_exponent = 4; % 路损系数 α
K = 0.01; % 路损常数 K
Power_UE = 24; % 用户的发射功率 24dbm    10^(24/10) mw;
Energy_loss_factor = 0.3; %电池的能量损失30%
T_slot = 1; % 一个时隙时常 1us   1秒=1000毫秒(ms)1秒=1000000 微秒(μs)
SINR_require = 8; % 需要的最低信噪比 10db


Num_DUE = 10; % DUE用户10个
 
OperaVol = 4; % 中继节点工作电压4V
% 可变参数
bandwidth = 720000; % 带宽720 000 hz = 720khz

 %% S和D之间距离对通信成功率的影响
 % D2D 通信坐标覆盖范围，假设发射方在0,300处，接收方在0,480处，距离180m
x_S = 0;
y_S = 280;
x_D = 0;

point_num = 10;
randomTotalTimeEnd = zeros(1,point_num);  % 随机方法总时间;
RateBaseTotalTimeEnd = zeros(1,point_num); % 以速率为基准的总时间
RateBufferBaseTotalTimeEnd = zeros(1,point_num); % 以速率、缓存空间为基准的总时间

randomTotalLossEnd = zeros(1,point_num); % 随机方法失败时间
RateBaseTotalLossEnd = zeros(1,point_num); % 以速率为基准失败时间
RateBufferBaseTotalLossEnd = zeros(1,point_num); % 以速率、缓存空间为基准的失败时间

for time =  1:1:1
    randomTotalTime = zeros(1,point_num);  % 随机方法总时间
    randomLoss = zeros(1,point_num);          % 随机方法失败时间

    RateBaseTotalTime = zeros(1,point_num); % 以速率为基准的总时间
    RateBaseLoss = zeros(1,point_num);          % 以速率为基准的失败时间

    RateBufferBaseTotalTime = zeros(1,point_num); % 以速率、缓存空间为基准的总时间
    RateBufferBaseTotalLoss = zeros(1,point_num); % 以速率、缓存空间为基准的失败时间

    local = 1;
    for y_D = 300:20:480 % 距离20-200 
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
        % 随机产生DUE 用户 10个
        x = 2*R_D2D*rand(1,10000) - R_D2D;
        y = 2*R_D2D*rand(1,10000) - R_D2D;
        x_DUE_temp = x(x.^2+y.^2<R_D2D*R_D2D);
        y_DUE_temp = y(x.^2+y.^2<R_D2D*R_D2D);
        x_DUE = x_DUE_temp(1:Num_DUE) + x_S;
        y_DUE = y_DUE_temp(1:Num_DUE) + y_S;
        scatter(x_DUE,y_DUE,'mh');
        hold on;
        axis([-500 500 -500 500]);


        % 生成中继节点的能量 暂且认为能量无限 0-2000 mAh 随机分布
        % 计算所有节点工作时间
        RelayEnery = 2000*rand(1,Num_DUE);
        t_residue = zeros(1,Num_DUE);
        for i = 1:1:Num_DUE
            t_residue(i) = 1000000 * judgeRelayWorkTime(RelayEnery(i),Power_UE,Energy_loss_factor,OperaVol); % us数
        end

        %方法一 随机选择链路 
        [tranTime,data_save_time,m1_random,m2_random] = randomSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        %方法二、三
        [RateBaseTime,RateBufferBaseTime,m1_RateBase,m2_RateBase,m1_RateBufferBase,m2_RateBufferBase] = comTwoWay(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);

        randomTotalTime(local) = tranTime;
        randomLoss(local) = m1_random + m2_random;

        RateBaseTotalTime(local) = RateBaseTime;
        RateBaseLoss(local) = m1_RateBase + m2_RateBase;

        RateBufferBaseTotalTime(local) = RateBufferBaseTime;
        RateBufferBaseTotalLoss(local) = m1_RateBufferBase + m2_RateBufferBase;
        local = local + 1
    end
    time
    randomTotalTimeEnd = randomTotalTimeEnd + randomTotalTime;
    RateBaseTotalTimeEnd = RateBaseTotalTimeEnd + RateBaseTotalTime;
    RateBufferBaseTotalTimeEnd = RateBufferBaseTotalTimeEnd + RateBufferBaseTotalTime;
    
    randomTotalLossEnd = randomTotalLossEnd + randomLoss;
    RateBaseTotalLossEnd = RateBaseTotalLossEnd + RateBaseLoss;
    RateBufferBaseTotalLossEnd = RateBufferBaseTotalLossEnd + RateBufferBaseTotalLoss;
end
% % 计算通信成功率
% randomSucc = (randomTotalTime - randomLoss)./randomTotalTime;
% RateBaseSucc = (RateBaseTotalTime - RateBaseLoss)./RateBaseTotalTime;
% RateBufferBaseSucc = (RateBufferBaseTotalTime - RateBufferBaseTotalLoss)./RateBufferBaseTotalTime;
randomTotalTimeEnd = randomTotalTimeEnd/time;
RateBaseTotalTimeEnd = RateBaseTotalTimeEnd/time;
RateBufferBaseTotalTimeEnd = RateBufferBaseTotalTimeEnd/time;

randomTotalLossEnd = randomTotalLossEnd/time;
RateBaseTotalLossEnd = RateBaseTotalLossEnd/time;
RateBufferBaseTotalLossEnd = RateBufferBaseTotalLossEnd/time;

ww = 20:20:200;
figure(2);
plot(ww,log(randomTotalTimeEnd),'cx--'); hold on;
plot(ww,log(RateBaseTotalTimeEnd),'bx--'); hold on;
plot(ww,log(RateBufferBaseTotalTimeEnd),'kx--'); hold on;

randomSucc = (randomTotalTimeEnd - randomTotalLossEnd)./randomTotalTimeEnd;
RateBaseSucc = (RateBaseTotalTimeEnd - RateBaseTotalLossEnd)./RateBaseTotalTimeEnd;
RateBufferBaseSucc = (RateBufferBaseTotalTimeEnd - RateBufferBaseTotalLossEnd)./RateBufferBaseTotalTimeEnd;

figure(3);
plot(ww,randomSucc,'cx--'); hold on;
plot(ww,RateBaseSucc,'bx--'); hold on;
plot(ww,RateBufferBaseSucc,'kx--'); hold on;
 