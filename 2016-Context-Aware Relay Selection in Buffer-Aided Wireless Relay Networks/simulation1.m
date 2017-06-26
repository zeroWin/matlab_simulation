clc,clear;
% 对原始论文的仿真

%% 基础变量
Packet_num = 50;  % 发送的包的个数
Relay_num = 4;    % 中继的个数
P_max = 10;         % mw 发送能够达到的最大功率
L_buffer = 10;      % 中继能够存储的个数
Noise = 1;            % mw 噪声功率
SNR = 10;            % 要满足的信噪比

%% 仿真1  固定中继节点的能量
% 变化发送节点的能量，计算两个图
% 图1：能够发送多少个包
% 图2：发送完50个包所需要的时间
E_relay = 600;      % mj 中继节点能量

% 仿能够发送多少个包

    
for E_source = 100:50:800 % 能量从100到800
    total_packet = 0;
    total_time = 0;
    for simulationTrials = 1:1:1000 % 仿真1000次
        Packect_s = Packet_num;
        Packect_d = 0;
        send_error = 0;
        Relay_contain = zeros(1,Relay_num);  % 初始化中继存储有多少个数据
        E_relay_con = zeros(1,Relay_num) +   E_relay; %  初始化中继剩余能量
        E_source_con = E_source;  % 发送节点剩余能量


        for testTime = 1:1:1000 %尝试选择1000次链路

            % 生成该时候的信道系数
            h = normrnd(0,1,1,Relay_num*2);

            % 计算s->r 的权重
            h_sr = h(1:Relay_num).^2; % s->r的信道系数
            P_sm = SNR*Noise./h_sr;     % 每个发送到每个中继所需要的最小功率
            W_sr = zeros(1,Relay_num);  % 保存每一条线路的权重
            e_m =    L_buffer -  Relay_contain; % 保存每个中继在t时刻还能存多少个节点
            s_m =    floor(E_relay_con ./ P_max) -Relay_contain;% 保存每个中继在t时刻还能发送多少个节点

            for i = 1:1:Relay_num % 计算每一条链路的权重
                if P_sm(i) <= min(P_max, E_source_con)  &&  Packect_s > 0
                    W_sr(i) =  E_source_con/P_sm(i) * min(e_m(i),s_m(i));  
                end
            end


            % 计算r->d 的权重
            h_rd = h(Relay_num+1:end).^2; % r->d的信道系数
            P_rd = SNR*Noise./h_rd;     % 每个中继发送所需要的最小功率 
            W_rd = zeros(1,Relay_num);  % 保存每一条线路的权重

            for i = 1:1:Relay_num % 计算每一条链路的权重
                if P_rd(i) <= min(P_max,E_relay_con(i))
                    W_rd(i) = E_relay_con(i) / P_rd(i) * Relay_contain(i);
                end
            end

            [maxNum, maxLocation] = max([W_sr,W_rd]);
            if(maxNum == 0) % 最大的权重为0 ，没有满足条件的链路 连续50次没有满足条件的链路，认为能量不足，无法完成发送
                send_error = send_error + 1;
            else
                if(maxLocation <= Relay_num) % 最大值在前面的链路
                     Relay_contain(maxLocation) = Relay_contain(maxLocation) + 1;
                     E_source_con  = E_source_con - P_sm(maxLocation); % 修改source剩余能量
                     Packect_s = Packect_s - 1;
                else % 最大值在后面的链路
                     maxLocation = maxLocation - Relay_num;
                     Relay_contain(maxLocation) = Relay_contain(maxLocation) - 1;
                     E_relay_con(maxLocation) = E_relay_con(maxLocation) -  P_rd(maxLocation);
                     Packect_d = Packect_d + 1;
                end        
                send_error = 0;
            end
            if(send_error == 50)
                %s = sprintf('发送无法完成，最后发送包成功接收个数为:%d\n,循环次数为:%d\n',Packect_d,testTime); 
                %disp(s);
                break;
            end
            if(Packect_d == Packet_num)
               % s = sprintf('发送成功循环次数为:%d\n',testTime); 
               % disp(s);
               total_time = total_time + testTime;
                break;
            end

        end
        total_packet  =   total_packet +  Packect_d;
        if( mod(simulationTrials,100) == 0)
            simulationTrials
        end
    end
    s = sprintf('平均接收个数为:%d\n',floor(total_packet/1000)); 
     disp(s);
    s = sprintf('平均接收耗费时间为:%d\n',floor(total_time/1000)); 
     disp(s);
end

    




