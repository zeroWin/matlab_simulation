% ������ �������ʺͻ���ѡ����·
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
function [tranTime,data_save_time] = RateBufferBaseSelect(t_residue,Relay_buffer,M_packet,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE)

% ����
SINR_require = 10;
R_min = bandwidth*log2(1 + 10);

data_save = [0 0 0 0 0 0 0 0 0 0]; % ��¼�ڵ��еĴ洢�������ĸ���
total_send = M_packet*M_packet_length;
total_receive = 0;

m1 = 0;
m2 = 0;
data_save_temp = zeros(1,10);
flag = 1;
for i = 1:1:1000000000
    array_fastFading =  exprnd(1,1,100000); % ����ָ���ֲ�����,�������ݣ�ȡ20��
    array_slowFading =  lognrnd(0,8,1,100000); % ������̬�����ֲ����飬�������ݣ�ȡ20��
    
    % ��������S�������м̶�CUE�û���SINR
    SINR_SToCue = judgeSINR_StoCue(x_S,y_S,x_CUE,y_CUE,array_fastFading(1:21),array_slowFading(1:21),bandwidth,Power_UE);
    SINR_RToCue = zeros(10,20);
    for relayNum = 1:1:10
        SINR_RToCue(relayNum,:) = judgeSINR_StoCue(x_DUE(relayNum),y_DUE(relayNum),x_CUE,y_CUE,array_fastFading(21*relayNum+1:21*relayNum+21),array_slowFading(21*relayNum+1:21*relayNum+21),bandwidth,Power_UE);
    end
    
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
    
    % ͬʱ����ÿ���ڵ��ʣ��������ʣ��ռ�ѡ�����ʺϵ���·����
    RateStoR = zeros(1,10);
    RateRtoD = zeros(1,10);
    % ����ÿ����·��Ȩ��
    W_StoR = zeros(1,10);
    W_RtoD = zeros(1,10);
    for relayNum = 1:1:10
        % S-Relay��·����
        if total_send > 0 % S�л���������
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
                RelaycanSend = R_min*t_residue(relayNum);
                % �ýڵ�ʣ��ռ仹�ܴ洢��������
                RelaycanSave = Relay_buffer * M_packet_length - data_save(relayNum);
                % ����Ȩ��
                W_StoR(relayNum) = min(RelaycanSend,RelaycanSave)*R;
                %W_StoR(relayNum) = R;
                % ��������
                RateStoR(relayNum) = R;
            end
        end
        % Relay->D��·����
        if data_save(relayNum) ~= 0 % �ڵ��д洢�������� ���Ż����ýڵ�
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
                W_RtoD(relayNum) = data_save(relayNum)*R;
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
            if total_send <= tranSize % ���ͽڵ����ݴ�����
                t_tran = total_send/R_min;
                 if data_save(maxLoc) + total_send <= Relay_buffer * M_packet_length && t_residue(maxLoc) >= t_tran*1000000
                     data_save(maxLoc) = data_save(maxLoc) + total_send;
                     total_send = 0;
                 else
                     m2 = m2 + 1;
                 end
            else
                t_tran = tranSize/R_min;
                if data_save(maxLoc) + tranSize <= Relay_buffer * M_packet_length && t_residue(maxLoc) >= t_tran*1000000
                    % ����Ӧ�м����������������Ҹ����м̽ڵ�洢��������
                    data_save(maxLoc) = data_save(maxLoc) + tranSize;
                    total_send = total_send - tranSize;
                else
                    m2 = m2 + 1;
                end
            end                        
        else % Relay->D��·
            maxLoc = maxLoc - 10;
            % ����һ��ʱ϶�����������
            tranSize = floor(RateRtoD(maxLoc) * (1/1000000));
            if tranSize >= data_save(maxLoc)  % ���Դ�������������ڴ洢��������
                % ���㴫����Ҫ���ѵ�ʱ��
                t_tran = data_save(maxLoc) / RateRtoD(maxLoc);
                %��������
                total_receive = total_receive + data_save(maxLoc); % �����ܽ��յ�������
                data_save(maxLoc) = 0;    % ���½ڵ��д洢��������
                t_residue(maxLoc) =  t_residue(maxLoc) - t_tran*1000000; % ���½ڵ�ʣ�๤��ʱ��           
            else
                 total_receive = total_receive + tranSize; % �����ܽ��յ�������
                 data_save(maxLoc) = data_save(maxLoc) - tranSize;    % ���½ڵ��д洢��������
                 t_residue(maxLoc) =  t_residue(maxLoc) - 1; % ���½ڵ�ʣ�๤��ʱ��       
            end                
        end
    else
        m1 = m1 + 1;
    end
    data_save_temp = [data_save_temp;data_save];
   if total_receive / 2048 > flag
       total_receive
       flag = flag + 2;
   end
    if total_receive >= M_packet*M_packet_length
        total_receive
        i
        m1
        m2
        break;
    end
    tranTime = i;
    data_save_time = data_save_temp([2:end],:);
end