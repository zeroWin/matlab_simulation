% �ȽϿ������ʵķ����� ͬʱ�������ʻ��淽���Ĵ���ʱ��
% ���룺
% 1. t_residue �ڵ�ʣ������
% 2.Relay_buffer �м��ܴ洢�İ��ĸ���
% 3.M_packet Ҫ����İ��ĸ���
% 4.M_packet_length ÿ������bit��
% 5.x_S,y_S ���ͽڵ�����
% 6.x_D,y_D ���սڵ�����
% 7.x_CUE,y_CUE ����CUE������
% 8.x_DUE,y_DUE ����DUE������
% 9.,bandwidth  ����
% 10.Power_UE ���书��

function [RateBaseTime,RateBufferBaseTime,totalDataRateBase,totalDataRateBufferBase] = comTwoWayTotalData(t_residue,Relay_buffer,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require)

    % ����
    R_min = bandwidth*log2(1 + 10^(SINR_require/10));

    % RateBase �õĲ���
    data_save_RateBase = [0 0 0 0 0 0 0 0 0 0]; % ��¼�ڵ��еĴ洢�������ĸ���
    total_send_RateBase = 1;
    total_receive_RateBase = 0;

    m1_RateBase = 0;
    m2_RateBase = 0;
    flag_RateBase = 1;
    
    % Rate Buffer Base �õĲ���
    data_save_RateBufferBase = [0 0 0 0 0 0 0 0 0 0]; % ��¼�ڵ��еĴ洢�������ĸ���
    total_send_RateBufferBase = 1;
    total_receive_RateBufferBase = 0;

    m1_RateBufferBase = 0;
    m2_RateBufferBase = 0;
    flag_RateBufferBase = 1;
    
    t_RateBase = t_residue;
    t_RateBufferBase = t_residue;
   
    for i = 1:1:1000000000000000000
        array_fastFading =  exprnd(1,1,2500); % ����ָ���ֲ�����,�������ݣ�ȡ20��
        array_slowFading =  lognrnd(0,8,1,2500); % ������̬�����ֲ����飬�������ݣ�ȡ20��
        
        % ��������S�������м̶�CUE�û���SINR
        SINR_SToCue = judgeSINR_StoCue(x_S,y_S,x_CUE,y_CUE,array_fastFading(1:21),array_slowFading(1:21),bandwidth,Power_UE);
        SINR_RToCue = zeros(10,20);
        for relayNum = 1:1:10
            SINR_RToCue(relayNum,:) = judgeSINR_StoCue(x_DUE(relayNum),y_DUE(relayNum),x_CUE,y_CUE,array_fastFading(21*relayNum+1:21*relayNum+21),array_slowFading(21*relayNum+1:21*relayNum+21),bandwidth,Power_UE);
        end

        % ��������CUE�û�SINR
        % S->Relay��SINR
        SINR_SToRelay= zeros(10,20);
        for relayNum = 1:1:10
            SINR_SToRelay(relayNum,:) = judgeSINR_CuetoD(x_S,y_S,x_DUE(relayNum),y_DUE(relayNum),x_CUE,y_CUE,array_fastFading(1000+21*(relayNum-1):1020+21*(relayNum-1)),array_slowFading(1000+21*(relayNum-1):1020+21*(relayNum-1)),bandwidth,Power_UE);
        end
       % Relay->D��SINR
       SINR_RelayToD = zeros(10,20);
        for relayNum = 1:1:10
            SINR_RelayToD(relayNum,:) = judgeSINR_CuetoD(x_DUE(relayNum),y_DUE(relayNum),x_D,y_D,x_CUE,y_CUE,array_fastFading([1500+relayNum-1,2000:2019]),array_slowFading([1500+relayNum-1,2000:2019]),bandwidth,Power_UE);
        end   

        %% ����2���������ʵķ���
        if (isempty(find(t_RateBase < 10))) % û�нڵ�����С��10
            % �Ҵ�������������·
            R_max = 0;
            R_maxloc = 0;
            for relayNum = 1:1:10
                % S-Relay��·����
                if total_send_RateBase > 0 % S�л���������
                    SINR_SToCueDB = 10*log10(SINR_SToCue);
                    SINR_SToRelayDB = 10*log10(SINR_SToRelay(relayNum,:));
                    % �ҵ���������ȵ��ŵ� 
                    canUsedChannel = intersect(find(SINR_SToCueDB > SINR_require),find(SINR_SToRelayDB > SINR_require));
                    SINR_temp = SINR_SToRelay(relayNum,:);
                    max_SINR= max(SINR_temp(canUsedChannel));
                    if(~isempty(max_SINR))  % ���������������ŵ���ȡ��������
                        % ���㴫������
                        R = bandwidth * log2(1 + max_SINR); % bit/s
                        if R > R_max
                            R_max = R;
                            R_maxloc = relayNum;
                        end
                    end
                end
                % Relay->D��·����
                if data_save_RateBase(relayNum) ~= 0 % �ڵ��д洢�������� ���Ż����ýڵ�
                    SINR_RToCueDB = SINR_RToCue(relayNum,:);
                    SINR_RelayToDDB = SINR_RelayToD(relayNum,:);
                    % �ҵ���������ȵ��ŵ� 
                    canUsedChannel = intersect(find(SINR_RToCueDB > SINR_require),find(SINR_RelayToDDB > SINR_require));
                    SINR_temp = SINR_RelayToD(relayNum,:);
                    max_SINR= max(SINR_temp(canUsedChannel));
                    if(~isempty(max_SINR))  % ���������������ŵ���ȡ��������
                        % ���㴫������
                        R = bandwidth * log2(1 + max_SINR); % bit/s
                        if R > R_max
                            R_max = R;
                            R_maxloc = relayNum + 10;
                        end
                    end        
                end
            end

            if R_max ~= 0 % ��Ϊ0��˵��������·����ͨ��Ҫ��
                % ����һ��ʱ϶�����������
                tranSize = floor(R_max * (1/1000000));
                if R_maxloc <= 10  % ǰ����·
                    t_tran = tranSize/R_min;
                    if data_save_RateBase(R_maxloc) + tranSize <= Relay_buffer * M_packet_length && t_RateBase(R_maxloc) >= t_tran*1000000
                        % ����Ӧ�м����������������Ҹ����м̽ڵ�洢��������
                        data_save_RateBase(R_maxloc) = data_save_RateBase(R_maxloc) + tranSize;
                    else
                        m2_RateBase = m2_RateBase + 1;
                    end        
                else % ������·
                    R_maxloc = R_maxloc - 10;
                    if tranSize >= data_save_RateBase(R_maxloc)  % ���Դ�������������ڴ洢��������
                        % ���㴫����Ҫ���ѵ�ʱ��
                        t_tran = data_save_RateBase(R_maxloc) / R_max;
                        %��������
                        total_receive_RateBase = total_receive_RateBase + data_save_RateBase(R_maxloc); % �����ܽ��յ�������
                        data_save_RateBase(R_maxloc) = 0;    % ���½ڵ��д洢��������
                        t_RateBase(R_maxloc) =  t_RateBase(R_maxloc) - t_tran*1000000; % ���½ڵ�ʣ�๤��ʱ��           
                    else
                         total_receive_RateBase = total_receive_RateBase + tranSize; % �����ܽ��յ�������
                         data_save_RateBase(R_maxloc) = data_save_RateBase(R_maxloc) - tranSize;    % ���½ڵ��д洢��������
                         t_RateBase(R_maxloc) =  t_RateBase(R_maxloc) - 1; % ���½ڵ�ʣ�๤��ʱ��       
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
        
        %% ���������������ʺ�buffer�ķ���
        if (isempty(find(t_RateBufferBase < 10))) % û�нڵ�����С��10
            % ͬʱ����ÿ���ڵ��ʣ��������ʣ��ռ�ѡ�����ʺϵ���·����
            RateStoR = zeros(1,10);
            RateRtoD = zeros(1,10);
            % ����ÿ����·��Ȩ��
            W_StoR = zeros(1,10);
            W_RtoD = zeros(1,10);
            for relayNum = 1:1:10
                % S-Relay��·����
                if total_send_RateBufferBase > 0 % S�л���������
                    SINR_SToCueDB = 10*log10(SINR_SToCue);
                    SINR_SToRelayDB = 10*log10(SINR_SToRelay(relayNum,:));
                    % �ҵ���������ȵ��ŵ� 
                    canUsedChannel = intersect(find(SINR_SToCueDB > SINR_require),find(SINR_SToRelayDB > SINR_require));
                    SINR_temp = SINR_SToRelay(relayNum,:);
                    max_SINR= max(SINR_temp(canUsedChannel));
                    if(~isempty(max_SINR))  % ���������������ŵ���ȡ��������
                        % ���㴫������
                        R = bandwidth * log2(1 + max_SINR); % bit/s
                        % ��ýڵ�ʣ�������ܹ����͵�������
                        RelaycanSend = R_min*t_RateBufferBase(relayNum)*(1/1000000);
                        % �ýڵ�ʣ��ռ仹�ܴ洢��������
                        RelaycanSave = Relay_buffer * M_packet_length - data_save_RateBufferBase(relayNum);
                        % ����Ȩ��
                        W_StoR(relayNum) = min(RelaycanSend,RelaycanSave)*R;
                        %W_StoR(relayNum) = R;
                        % ��������
                        RateStoR(relayNum) = R;
                    end
                end
                % Relay->D��·����
                if data_save_RateBufferBase(relayNum) ~= 0 % �ڵ��д洢�������� ���Ż����ýڵ�
                    SINR_RToCueDB = SINR_RToCue(relayNum,:);
                    SINR_RelayToDDB = SINR_RelayToD(relayNum,:);
                    % �ҵ���������ȵ��ŵ� 
                    canUsedChannel = intersect(find(SINR_RToCueDB > SINR_require),find(SINR_RelayToDDB > SINR_require));
                    SINR_temp = SINR_RelayToD(relayNum,:);
                    max_SINR= max(SINR_temp(canUsedChannel));
                    if(~isempty(max_SINR))  % ���������������ŵ���ȡ��������
                        % ���㴫������
                        R = bandwidth * log2(1 + max_SINR); % bit/s

                        % ����Ȩ�� �Ѿ��洢��������*����
                        W_RtoD(relayNum) = data_save_RateBufferBase(relayNum)*R;
                        %W_RtoD(relayNum) = R;
                        % ��������
                        RateRtoD(relayNum) = R;                
                    end        
                end
            end

            % ѡ��Ȩ��������·���д���
            [maxWeight,maxLoc] = max([W_StoR W_RtoD]);
            if(maxWeight ~= 0) % ���ڿ��Է��͵Ľڵ�
                if maxLoc <= 10 % S->Relay ��·
                    % ����һ��ʱ϶�����������
                    tranSize = floor(RateStoR(maxLoc) * (1/1000000));
                    t_tran = tranSize/R_min;
                    if data_save_RateBufferBase(maxLoc) + tranSize <= Relay_buffer * M_packet_length && t_RateBufferBase(maxLoc) >= t_tran*1000000
                        % ����Ӧ�м����������������Ҹ����м̽ڵ�洢��������
                        data_save_RateBufferBase(maxLoc) = data_save_RateBufferBase(maxLoc) + tranSize;
                    else
                        m2_RateBufferBase = m2_RateBufferBase + 1;
                    end                       
                else % Relay->D��·
                    maxLoc = maxLoc - 10;
                    % ����һ��ʱ϶�����������
                    tranSize = floor(RateRtoD(maxLoc) * (1/1000000));
                    if tranSize >= data_save_RateBufferBase(maxLoc)  % ���Դ�������������ڴ洢��������
                        % ���㴫����Ҫ���ѵ�ʱ��
                        t_tran = data_save_RateBufferBase(maxLoc) / RateRtoD(maxLoc);
                        %��������
                        total_receive_RateBufferBase = total_receive_RateBufferBase + data_save_RateBufferBase(maxLoc); % �����ܽ��յ�������
                        data_save_RateBufferBase(maxLoc) = 0;    % ���½ڵ��д洢��������
                        t_RateBufferBase(maxLoc) =  t_RateBufferBase(maxLoc) - t_tran*1000000; % ���½ڵ�ʣ�๤��ʱ��           
                    else
                         total_receive_RateBufferBase = total_receive_RateBufferBase + tranSize; % �����ܽ��յ�������
                         data_save_RateBufferBase(maxLoc) = data_save_RateBufferBase(maxLoc) - tranSize;    % ���½ڵ��д洢��������
                         t_RateBufferBase(maxLoc) =  t_RateBufferBase(maxLoc) - 1; % ���½ڵ�ʣ�๤��ʱ��       
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