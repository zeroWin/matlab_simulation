clc,clear;
%% 场景图仿真
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
Relay_buffer = 5; % 能存储包个数5
Noise_density = -174; % 噪声功率谱密度 -174dbm/hz    dbm = 10*lg(mw)
Path_loss_exponent = 4; % 路损系数 α
K = 0.01; % 路损常数 K
Power_UE = 24; % 用户的发射功率 24dbm    10^(24/10) mw;
Energy_loss_factor = 0.3; %电池的能量损失30%
T_slot = 1; % 一个时隙时常 1us   1秒=1000毫秒(ms)1秒=1000000 微秒(μs)
SINR_require = 10; % 需要的最低信噪比 10db

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

%% 方法一 随机选择链路 
% array_fastFading =  exprnd(1);  % 均值为1的指数分布
% array_slowFading = lognrnd(0,8); %均值为0，方差为8db的正态对数分布

% [tranTime,data_save_time] = randomSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE);
% 
% %% 方法二 以速率为准则选择链路
% %[tranTime,data_save_time] = RateBaseSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE);
% 
% %% 方法三 同时考虑速率，能量，缓存空间选择链路
% %[tranTime,data_save_time] = RateBufferBaseSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE);
% 
% % 比较两种方法传输
%  [RateBaseTime,RateBufferBaseTime] = comTwoWay(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE);
% 
%  
 

 
 
 