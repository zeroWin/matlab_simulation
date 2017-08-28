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
SINR_require = 20; % 需要的最低信噪比 20db


Num_DUE = 10; % DUE用户10个
 
OperaVol = 4; % 中继节点工作电压4V
% 可变参数
bandwidth = 720000; % 带宽720 000 hz = 720khz
%bandwidth = 1; % 带宽720 000 hz = 720khz

 %% S和D之间距离对通信成功率的影响
 % D2D 通信坐标覆盖范围，假设发射方在0,250处，接收方在0,470处，距离220m
x_S = 0;
y_S = 250;
x_D = 0;

point_num = 6;

randomTotalTimeEnd = zeros(1,point_num);  % 随机方法总时间;
RateBaseTotalTimeEnd = zeros(1,point_num); % 以速率为基准的总时间
RateBufferBaseTotalTimeEnd = zeros(1,point_num); % 以速率、缓存空间为基准的总时间

randomTotalData = zeros(1,point_num);  % 随机方法总时间;
RateBaseTotalData = zeros(1,point_num); % 以速率为基准的总时间
RateBufferBaseTotalData = zeros(1,point_num); % 以速率、缓存空间为基准的总时间
t_save = zeros(1,Num_DUE);
for time =  1:1:20
    
    local = 1;
    for y_D = 270:40:470 % 距离20-220  270-40-470
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
        %RelayEnery = 0.01*rand(1,Num_DUE);
        RelayEnery = 0.0001*ones(1,10);
        t_residue = zeros(1,Num_DUE);
        for i = 1:1:Num_DUE
            t_residue(i) = 1000000*judgeRelayWorkTime(RelayEnery(i),Power_UE,Energy_loss_factor,OperaVol); % us数
        end
        
        %方法一 随机选择链路 
        [tranTime,totalData] = randomSelectTotalData(t_residue,Relay_buffer,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);
        %方法二、三
        [RateBaseTime,RateBufferBaseTime,totalDataRateBase,totalDataRateBufferBase] = comTwoWayTotalData(t_residue,Relay_buffer,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require);

        randomTotalTimeEnd(local) =  randomTotalTimeEnd(local) + tranTime;
        RateBaseTotalTimeEnd(local) = RateBaseTotalTimeEnd(local) +  RateBaseTime;
        RateBufferBaseTotalTimeEnd(local) =  RateBufferBaseTotalTimeEnd(local) + RateBufferBaseTime;

        randomTotalData(local) = randomTotalData(local) + totalData;
        RateBaseTotalData(local) = RateBaseTotalData(local) +totalDataRateBase;
        RateBufferBaseTotalData(local) =  RateBufferBaseTotalData(local) + totalDataRateBufferBase;
        
        local = local + 1
    end
    time
end
% % 计算通信成功率
% randomSucc = (randomTotalTime - randomLoss)./randomTotalTime;
% RateBaseSucc = (RateBaseTotalTime - RateBaseLoss)./RateBaseTotalTime;
% RateBufferBaseSucc = (RateBufferBaseTotalTime - RateBufferBaseTotalLoss)./RateBufferBaseTotalTime;
% 8db
randomTotalTimeEnd = randomTotalTimeEnd/time;
RateBaseTotalTimeEnd = RateBaseTotalTimeEnd/time;
RateBufferBaseTotalTimeEnd = RateBufferBaseTotalTimeEnd/time;

randomTotalData = randomTotalData/time;
RateBaseTotalData = RateBaseTotalData/time;
RateBufferBaseTotalData = RateBufferBaseTotalData/time;

% 画图
ww = 20:40:220;
%ww = 20:20:80;
figure(2);
plot(ww,log(randomTotalTimeEnd),'c^--'); hold on;
plot(ww,log(RateBaseTotalTimeEnd),'mo--'); hold on;
plot(ww,log(RateBufferBaseTotalTimeEnd),'rs--'); hold on;

figure(3);
plot(ww,(randomTotalData),'c^--'); hold on;
plot(ww,(RateBaseTotalData),'mo--'); hold on;
plot(ww,(RateBufferBaseTotalData),'rs--'); hold on;

plot(ww,(randomTotalTimeEnd10db),'b^-'); hold on;
plot(ww,(RateBaseTotalTimeEnd10db),'ko-'); hold on;
plot(ww,(RateBufferBaseTotalTimeEnd10db),'gs-'); hold on;
