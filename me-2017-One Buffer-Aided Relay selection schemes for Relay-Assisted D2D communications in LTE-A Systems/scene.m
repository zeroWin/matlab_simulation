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
Relay_buffer = 10; % 能存储包个数10
Noise_density = -174; % 噪声功率谱密度 -174dbm/hz    dbm = 10*lg(mw)
Path_loss_exponent = 4; % 路损系数 α
K = 0.01; % 路损常数 K
Power_UE = 24; % 用户的发射功率 24dbm    10^(24/10) mw;
Energy_loss_factor = 0.3; %电池的能量损失30%
T_slot = 1; % 一个时隙时常 1us   1秒=1000毫秒(ms)1秒=1000000 微秒(μs)
SINR_require = 10; % 需要的最低信噪比 10db

Num_CUE = 20; % CUE用户20个
Num_DUE = 10; % DUE用户10个
 

% 可变参数
bandwidth = 720000; % 带宽720 000 hz = 720khz
R_min = bandwidth*log2(1 + 10);
% 生成中继节点的能量 暂且认为能量无限


%% 方法一 随机选择链路 
% array_fastFading =  exprnd(1);  % 均值为1的指数分布
% array_slowFading = lognrnd(0,8); %均值为0，方差为8db的正态对数分布


road_array = [1 2 3 4 5 6 7 8 9 10 0 0 0 0 0 0 0 0 0 0]; % 记录那些链路能选择
data_save = [0 0 0 0 0 0 0 0 0 0]; % 记录节点中的存储数据量的个数
road_array_temp = find(road_array > 0);
index = randperm(length(road_array_temp));
road_select = road_array_temp(index(1)); % 选择出的路径

% 计算该链路的信噪比
% 1. 选择信噪比最高的那个复用信道
% 2. 计算传输速率，是否满足最低信噪比，满足进行下一步，不满足下一个时隙重新选择
% 3. 计算传输的数据量
% 4. 根据选择的是前向链路还是后向链路分情况处理
% 情况1：如果选择路径是前面的链路，计算节点能量以最慢速率是否足够传输？计算剩余空间是否足够？
% 足够，就传，不够，下一个时隙重新选择链路
% 如果选择的是后向链路，传输选择节点中的数据。

array_fastFading =  exprnd(1,1,100000); % 生成指数分布数组,大量数据，取20个
array_slowFading =  lognrnd(0,8,1,100000); % 生成正态对数分布数组，大量数据，取20个

% 获取选择路径的DUE的坐标
if road_select > 10  % Relay->D
    Relay_x = x_DUE(road_select-10);
    Relay_y = y_DUE(road_select-10);
 
    
else  % S -> Relay
    Relay_x = x_DUE(road_select);
    Relay_y = y_DUE(road_select); 
    
    % 求每个复用信道下CUE用户对选择中继节点的干扰
    [SINR_CueToRelay, SINR_RelayToCue] = judgeSINR_StoRelay(Relay_x,Relay_y,x_S,y_S,x_CUE,y_CUE,array_fastFading,array_slowFading,bandwidth,Power_UE);
    SINR_CueToRelayDB = 10*log10(SINR_CueToRelay);
    SINR_RelayToCueDB = 10*log10(SINR_RelayToCue);

    % 找到满足信噪比的信道 
    canUsedChannel = intersect(find(SINR_CueToRelayDB > SINR_require),find(SINR_RelayToCueDB> SINR_require));
    max_SINR= max(SINR_CueToRelay(canUsedChannel));
    if(length(max_SINR))  % 存在满足条件的信道，取最大信噪比
        % 计算传输速率
        R = bandwidth * log2(1 + max_SINR); % bit/s
        % 计算一个时隙能传输的数据量
        tranSize = R * (1/1000000);
        
        
        
        % 判断中继节点能量和空间是否足够
        % 计算节点剩余能量可以工作的时间t_residue
        % 计算节点以最慢速率传输这些数据量需要的时间t_tran s
        t_residue = 1111;
        t_tran = tranSize/R_min;
        if data_save(road_select) + tranSize <= Relay_buffer * M_packet && t_residue >= t_tran
            % 给相应中继增加数据量，并且更新中继节点存储的数据量
            data_save(road_select) = data_save(road_select) + tranSize;
            road_array(road_select + 10) = road_select + 10;
        end
    end
end



 %% 方法二 以速率为准则选择链路


%% 方法三 同时考虑速率，能量，缓存空间选择链路



