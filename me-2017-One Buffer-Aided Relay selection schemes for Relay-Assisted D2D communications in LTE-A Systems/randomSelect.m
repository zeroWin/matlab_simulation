% 方法一 随机选择链路 
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
function [tranTime,data_save_time] = randomSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE)

% 常数
SINR_require = 10;
R_min = bandwidth*log2(1 + 10);

road_array = [1 2 3 4 5 6 7 8 9 10 0 0 0 0 0 0 0 0 0 0]; % 记录那些链路能选择
data_save = [0 0 0 0 0 0 0 0 0 0]; % 记录节点中的存储数据量的个数
total_send = M_packet*M_packet_length;
total_receive = 0;

m1 = 0;
m2 = 0;
data_save_temp = zeros(1,10);
flag = 1;
for i = 1:1:1000000000
    road_array_temp = find(road_array > 0);
    index = randperm(length(road_array_temp));
    road_select = road_array_temp(index(1)); % 选择出的路径


    % 需要多少个时隙可以发完所有的数据
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


        % 求每个复用信道下CUE用户对选择中继节点的干扰
        [SINR_CueToD, SINR_RelayToCue] = judgeSINR_RelaytoD(Relay_x,Relay_y,x_D,y_D,x_CUE,y_CUE,array_fastFading,array_slowFading,bandwidth,Power_UE);
        SINR_CueToD = 10*log10(SINR_CueToD);
        SINR_RelayToCue = 10*log10(SINR_RelayToCue); 

        % 找到满足信噪比的信道
        canUsedChannel = intersect(find(SINR_CueToD > SINR_require),find(SINR_RelayToCue> SINR_require));
        max_SINR= max(SINR_CueToD(canUsedChannel));

        if(~isempty(max_SINR))  % 存在满足条件的信道，取最大信噪比
            % 计算传输速率
            R = bandwidth * log2(1 + max_SINR); % bit/s
            % 计算一个时隙能传输的数据量
            tranSize = floor(R * (1/1000000));

            if tranSize >= data_save(road_select-10)  % 可以传输的数据量大于存储的数据量
                % 计算传输需要花费的时间
                t_tran = data_save(road_select-10) / R;
                %更新数据
                total_receive = total_receive + data_save(road_select-10); % 更新总接收到的数据
                data_save(road_select - 10) = 0;    % 更新节点中存储的数据量
                road_array(road_select) = 0;
                t_residue(road_select - 10) =  t_residue(road_select - 10) - t_tran*1000000; % 更新节点剩余工作时间           
            else
                 total_receive = total_receive + tranSize; % 更新总接收到的数据
                 data_save(road_select - 10) = data_save(road_select - 10) - tranSize;    % 更新节点中存储的数据量
                 t_residue(road_select - 10) =  t_residue(road_select - 10) - 1; % 更新节点剩余工作时间       
            end
        else
            m1 = m1 + 1;
        end
    else  % S -> Relay
        Relay_x = x_DUE(road_select);
        Relay_y = y_DUE(road_select); 

        % 求每个复用信道下CUE用户对选择中继节点的干扰
        [SINR_CueToRelay, SINR_SToCue] = judgeSINR_StoRelay(Relay_x,Relay_y,x_S,y_S,x_CUE,y_CUE,array_fastFading,array_slowFading,bandwidth,Power_UE);
        SINR_CueToRelayDB = 10*log10(SINR_CueToRelay);
        SINR_SToCueDB = 10*log10(SINR_SToCue);

        % 找到满足信噪比的信道 
        canUsedChannel = intersect(find(SINR_CueToRelayDB > SINR_require),find(SINR_SToCueDB> SINR_require));
        max_SINR= max(SINR_CueToRelay(canUsedChannel));
        if(~isempty(max_SINR))  % 存在满足条件的信道，取最大信噪比
            % 计算传输速率
            R = bandwidth * log2(1 + max_SINR); % bit/s
            % 计算一个时隙能传输的数据量
            tranSize = floor(R * (1/1000000));

            % 判断中继节点能量和空间是否足够
            % 计算节点剩余能量可以工作的时间t_residue s
            % 计算节点以最慢速率传输这些数据量需要的时间t_tran s
            if total_send <= tranSize % 发送节点数据传完了
                t_tran = total_send/R_min;
                 if data_save(road_select) + total_send <= Relay_buffer * M_packet_length && t_residue(road_select) >= t_tran*1000000
                     road_array(1:10) = 0;
                     data_save(road_select) = data_save(road_select) + total_send;
                     road_array(road_select + 10) = road_select + 10;
                     total_send = 0;
                 else
                     m2 = m2 + 1;
                 end
            else
                t_tran = tranSize/R_min;
                if data_save(road_select) + tranSize <= Relay_buffer * M_packet_length && t_residue(road_select) >= t_tran*1000000
                    % 给相应中继增加数据量，并且更新中继节点存储的数据量
                    data_save(road_select) = data_save(road_select) + tranSize;
                    road_array(road_select + 10) = road_select + 10;
                    total_send = total_send - tranSize;
                 else
                     m2 = m2 + 1;                    
                end
            end
        else
            m1 = m1 + 1;
        end
    end
   data_save_temp = [data_save_temp;data_save];
   if total_receive / 1024 > flag
       total_receive
       flag = flag + 1;
   end
    if total_receive >= M_packet*M_packet_length
        total_receive
        i
        m1
        m2
        break;
    end
end
    tranTime = i;
    data_save_time = data_save_temp([2:end],:);
end