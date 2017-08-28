% 比较考虑速率的方法和 同时考虑速率缓存方法的传输时间
% 输入：
% 1. t_residue 节点剩余能量
% 2.Relay_buffer 中继能存储的包的个数
% 3.M_packet 要传输的包的个数
% 4.M_packet_length 每个包的bit数
% 5.x_S,y_S 发送节点坐标
% 6.x_D,y_D 接收节点坐标
% 7.x_CUE,y_CUE 所有CUE的坐标
% 8.x_DUE,y_DUE 所有DUE的坐标
% 9.,bandwidth  带宽
% 10.Power_UE 发射功率

function [RateBaseTime,RateBufferBaseTime,totalDataRateBase,totalDataRateBufferBase] = comTwoWayTotalData(t_residue,Relay_buffer,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require)

    % 常数
    R_min = bandwidth*log2(1 + 10^(SINR_require/10));

    % RateBase 用的参数
    data_save_RateBase = [0 0 0 0 0 0 0 0 0 0]; % 记录节点中的存储数据量的个数
    total_send_RateBase = 1;
    total_receive_RateBase = 0;

    m1_RateBase = 0;
    m2_RateBase = 0;
    flag_RateBase = 1;
    
    % Rate Buffer Base 用的参数
    data_save_RateBufferBase = [0 0 0 0 0 0 0 0 0 0]; % 记录节点中的存储数据量的个数
    total_send_RateBufferBase = 1;
    total_receive_RateBufferBase = 0;

    m1_RateBufferBase = 0;
    m2_RateBufferBase = 0;
    flag_RateBufferBase = 1;
    
    t_RateBase = t_residue;
    t_RateBufferBase = t_residue;
   
    for i = 1:1:1000000000000000000
        array_fastFading =  exprnd(1,1,2500); % 生成指数分布数组,大量数据，取20个
        array_slowFading =  lognrnd(0,8,1,2500); % 生成正态对数分布数组，大量数据，取20个
        
        % 计算所有S和所有中继对CUE用户的SINR
        SINR_SToCue = judgeSINR_StoCue(x_S,y_S,x_CUE,y_CUE,array_fastFading(1:21),array_slowFading(1:21),bandwidth,Power_UE);
        SINR_RToCue = zeros(10,20);
        for relayNum = 1:1:10
            SINR_RToCue(relayNum,:) = judgeSINR_StoCue(x_DUE(relayNum),y_DUE(relayNum),x_CUE,y_CUE,array_fastFading(21*relayNum+1:21*relayNum+21),array_slowFading(21*relayNum+1:21*relayNum+21),bandwidth,Power_UE);
        end

        % 计算所有CUE用户SINR
        % S->Relay的SINR
        SINR_SToRelay= zeros(10,20);
        for relayNum = 1:1:10
            SINR_SToRelay(relayNum,:) = judgeSINR_CuetoD(x_S,y_S,x_DUE(relayNum),y_DUE(relayNum),x_CUE,y_CUE,array_fastFading(1000+21*(relayNum-1):1020+21*(relayNum-1)),array_slowFading(1000+21*(relayNum-1):1020+21*(relayNum-1)),bandwidth,Power_UE);
        end
       % Relay->D的SINR
       SINR_RelayToD = zeros(10,20);
        for relayNum = 1:1:10
            SINR_RelayToD(relayNum,:) = judgeSINR_CuetoD(x_DUE(relayNum),y_DUE(relayNum),x_D,y_D,x_CUE,y_CUE,array_fastFading([1500+relayNum-1,2000:2019]),array_slowFading([1500+relayNum-1,2000:2019]),bandwidth,Power_UE);
        end   

        %% 方法2：基于速率的方法
        if (isempty(find(t_RateBase < 10))) % 没有节点能量小于10
            % 找传输速率最大的链路
            R_max = 0;
            R_maxloc = 0;
            for relayNum = 1:1:10
                % S-Relay链路速率
                if total_send_RateBase > 0 % S中还存在数据
                    SINR_SToCueDB = 10*log10(SINR_SToCue);
                    SINR_SToRelayDB = 10*log10(SINR_SToRelay(relayNum,:));
                    % 找到满足信噪比的信道 
                    canUsedChannel = intersect(find(SINR_SToCueDB > SINR_require),find(SINR_SToRelayDB > SINR_require));
                    SINR_temp = SINR_SToRelay(relayNum,:);
                    max_SINR= max(SINR_temp(canUsedChannel));
                    if(~isempty(max_SINR))  % 存在满足条件的信道，取最大信噪比
                        % 计算传输速率
                        R = bandwidth * log2(1 + max_SINR); % bit/s
                        if R > R_max
                            R_max = R;
                            R_maxloc = relayNum;
                        end
                    end
                end
                % Relay->D链路速率
                if data_save_RateBase(relayNum) ~= 0 % 节点中存储的有数据 ，才会计算该节点
                    SINR_RToCueDB = SINR_RToCue(relayNum,:);
                    SINR_RelayToDDB = SINR_RelayToD(relayNum,:);
                    % 找到满足信噪比的信道 
                    canUsedChannel = intersect(find(SINR_RToCueDB > SINR_require),find(SINR_RelayToDDB > SINR_require));
                    SINR_temp = SINR_RelayToD(relayNum,:);
                    max_SINR= max(SINR_temp(canUsedChannel));
                    if(~isempty(max_SINR))  % 存在满足条件的信道，取最大信噪比
                        % 计算传输速率
                        R = bandwidth * log2(1 + max_SINR); % bit/s
                        if R > R_max
                            R_max = R;
                            R_maxloc = relayNum + 10;
                        end
                    end        
                end
            end

            if R_max ~= 0 % 不为0，说明存在链路满足通信要求
                % 计算一个时隙传输的数据数
                tranSize = floor(R_max * (1/1000000));
                if R_maxloc <= 10  % 前向链路
                    t_tran = tranSize/R_min;
                    if data_save_RateBase(R_maxloc) + tranSize <= Relay_buffer * M_packet_length && t_RateBase(R_maxloc) >= t_tran*1000000
                        % 给相应中继增加数据量，并且更新中继节点存储的数据量
                        data_save_RateBase(R_maxloc) = data_save_RateBase(R_maxloc) + tranSize;
                    else
                        m2_RateBase = m2_RateBase + 1;
                    end        
                else % 后向链路
                    R_maxloc = R_maxloc - 10;
                    if tranSize >= data_save_RateBase(R_maxloc)  % 可以传输的数据量大于存储的数据量
                        % 计算传输需要花费的时间
                        t_tran = data_save_RateBase(R_maxloc) / R_max;
                        %更新数据
                        total_receive_RateBase = total_receive_RateBase + data_save_RateBase(R_maxloc); % 更新总接收到的数据
                        data_save_RateBase(R_maxloc) = 0;    % 更新节点中存储的数据量
                        t_RateBase(R_maxloc) =  t_RateBase(R_maxloc) - t_tran*1000000; % 更新节点剩余工作时间           
                    else
                         total_receive_RateBase = total_receive_RateBase + tranSize; % 更新总接收到的数据
                         data_save_RateBase(R_maxloc) = data_save_RateBase(R_maxloc) - tranSize;    % 更新节点中存储的数据量
                         t_RateBase(R_maxloc) =  t_RateBase(R_maxloc) - 1; % 更新节点剩余工作时间       
                    end            
                end
            else
                m1_RateBase = m1_RateBase +1;
            end

           if total_receive_RateBase / 102400 > flag_RateBase
               total_receive_RateBase
               t_RateBase
               flag_RateBase = flag_RateBase + 1;
           end
        elseif total_send_RateBase == 1
            total_send_RateBase = 2;
            RateBaseTime = i;
            totalDataRateBase = total_receive_RateBase;
        end
        
        %% 方法三：基于速率和buffer的方法
        if (isempty(find(t_RateBufferBase < 10))) % 没有节点能量小于10
            % 同时考虑每个节点的剩余能量和剩余空间选出最适合的链路发送
            RateStoR = zeros(1,10);
            RateRtoD = zeros(1,10);
            % 计算每跳链路的权重
            W_StoR = zeros(1,10);
            W_RtoD = zeros(1,10);
            for relayNum = 1:1:10
                % S-Relay链路速率
                if total_send_RateBufferBase > 0 % S中还存在数据
                    SINR_SToCueDB = 10*log10(SINR_SToCue);
                    SINR_SToRelayDB = 10*log10(SINR_SToRelay(relayNum,:));
                    % 找到满足信噪比的信道 
                    canUsedChannel = intersect(find(SINR_SToCueDB > SINR_require),find(SINR_SToRelayDB > SINR_require));
                    SINR_temp = SINR_SToRelay(relayNum,:);
                    max_SINR= max(SINR_temp(canUsedChannel));
                    if(~isempty(max_SINR))  % 存在满足条件的信道，取最大信噪比
                        % 计算传输速率
                        R = bandwidth * log2(1 + max_SINR); % bit/s
                        % 求该节点剩余能量能够发送的数据量
                        RelaycanSend = R_min*t_RateBufferBase(relayNum)*(1/1000000);
                        % 该节点剩余空间还能存储的数据量
                        RelaycanSave = Relay_buffer * M_packet_length - data_save_RateBufferBase(relayNum);
                        % 计算权重
                        W_StoR(relayNum) = min(RelaycanSend,RelaycanSave)*R;
                        %W_StoR(relayNum) = R;
                        % 保存速率
                        RateStoR(relayNum) = R;
                    end
                end
                % Relay->D链路速率
                if data_save_RateBufferBase(relayNum) ~= 0 % 节点中存储的有数据 ，才会计算该节点
                    SINR_RToCueDB = SINR_RToCue(relayNum,:);
                    SINR_RelayToDDB = SINR_RelayToD(relayNum,:);
                    % 找到满足信噪比的信道 
                    canUsedChannel = intersect(find(SINR_RToCueDB > SINR_require),find(SINR_RelayToDDB > SINR_require));
                    SINR_temp = SINR_RelayToD(relayNum,:);
                    max_SINR= max(SINR_temp(canUsedChannel));
                    if(~isempty(max_SINR))  % 存在满足条件的信道，取最大信噪比
                        % 计算传输速率
                        R = bandwidth * log2(1 + max_SINR); % bit/s

                        % 计算权重 已经存储的数据数*速率
                        W_RtoD(relayNum) = data_save_RateBufferBase(relayNum)*R;
                        %W_RtoD(relayNum) = R;
                        % 保存速率
                        RateRtoD(relayNum) = R;                
                    end        
                end
            end

            % 选择权重最大的链路进行传输
            [maxWeight,maxLoc] = max([W_StoR W_RtoD]);
            if(maxWeight ~= 0) % 存在可以发送的节点
                if maxLoc <= 10 % S->Relay 链路
                    % 计算一个时隙传输的数据数
                    tranSize = floor(RateStoR(maxLoc) * (1/1000000));
                    t_tran = tranSize/R_min;
                    if data_save_RateBufferBase(maxLoc) + tranSize <= Relay_buffer * M_packet_length && t_RateBufferBase(maxLoc) >= t_tran*1000000
                        % 给相应中继增加数据量，并且更新中继节点存储的数据量
                        data_save_RateBufferBase(maxLoc) = data_save_RateBufferBase(maxLoc) + tranSize;
                    else
                        m2_RateBufferBase = m2_RateBufferBase + 1;
                    end                       
                else % Relay->D链路
                    maxLoc = maxLoc - 10;
                    % 计算一个时隙传输的数据数
                    tranSize = floor(RateRtoD(maxLoc) * (1/1000000));
                    if tranSize >= data_save_RateBufferBase(maxLoc)  % 可以传输的数据量大于存储的数据量
                        % 计算传输需要花费的时间
                        t_tran = data_save_RateBufferBase(maxLoc) / RateRtoD(maxLoc);
                        %更新数据
                        total_receive_RateBufferBase = total_receive_RateBufferBase + data_save_RateBufferBase(maxLoc); % 更新总接收到的数据
                        data_save_RateBufferBase(maxLoc) = 0;    % 更新节点中存储的数据量
                        t_RateBufferBase(maxLoc) =  t_RateBufferBase(maxLoc) - t_tran*1000000; % 更新节点剩余工作时间           
                    else
                         total_receive_RateBufferBase = total_receive_RateBufferBase + tranSize; % 更新总接收到的数据
                         data_save_RateBufferBase(maxLoc) = data_save_RateBufferBase(maxLoc) - tranSize;    % 更新节点中存储的数据量
                         t_RateBufferBase(maxLoc) =  t_RateBufferBase(maxLoc) - 1; % 更新节点剩余工作时间       
                    end                
                end
            else
                m1_RateBufferBase = m1_RateBufferBase + 1;
            end

            if total_receive_RateBufferBase / 102400 > flag_RateBufferBase
               total_receive_RateBufferBase
               t_RateBufferBase
               flag_RateBufferBase = flag_RateBufferBase + 1;
            end
        elseif total_send_RateBufferBase == 1
            total_send_RateBufferBase = 2;
            RateBufferBaseTime = i;
            totalDataRateBufferBase = total_receive_RateBufferBase;
        end
        if total_send_RateBase == 2 && total_send_RateBufferBase == 2
            total_receive_RateBase
            t_RateBase
            total_receive_RateBufferBase
            t_RateBufferBase
           break;
        end
    end
end