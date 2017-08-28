% ����һ ���ѡ����· 
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
function [tranTime,totalData] = randomSelectTotalData(t_residue,Relay_buffer,M_packet_length,x_S,y_S,x_D,y_D,x_CUE,y_CUE,x_DUE,y_DUE,bandwidth,Power_UE,SINR_require)

% ����
R_min = bandwidth*log2(1 + 10^(SINR_require/10));

road_array = [1 2 3 4 5 6 7 8 9 10 0 0 0 0 0 0 0 0 0 0]; % ��¼��Щ��·��ѡ��
data_save = [0 0 0 0 0 0 0 0 0 0]; % ��¼�ڵ��еĴ洢�������ĸ���
total_receive = 0;

m1 = 0;
m2 = 0;

flag = 1;
for i = 1:1:1000000000000000000
    road_array_temp = find(road_array > 0);
    index = randperm(length(road_array_temp));
    road_select = road_array_temp(index(1)); % ѡ�����·��


    % ��Ҫ���ٸ�ʱ϶���Է������е�����
    % �������·�������
    % 1. ѡ���������ߵ��Ǹ������ŵ�
    % 2. ���㴫�����ʣ��Ƿ������������ȣ����������һ������������һ��ʱ϶����ѡ��
    % 3. ���㴫���������
    % 4. ����ѡ�����ǰ����·���Ǻ�����·���������
    % ���1�����ѡ��·����ǰ�����·������ڵ����������������Ƿ��㹻���䣿����ʣ��ռ��Ƿ��㹻��
    % �㹻���ʹ�����������һ��ʱ϶����ѡ����·
    % ���ѡ����Ǻ�����·������ѡ��ڵ��е����ݡ�

    array_fastFading =  exprnd(1,1,1000); % ����ָ���ֲ�����,�������ݣ�ȡ20��
    array_slowFading =  lognrnd(0,8,1,1000); % ������̬�����ֲ����飬�������ݣ�ȡ20��

    % ��ȡѡ��·����DUE������
    if road_select > 10  % Relay->D
        Relay_x = x_DUE(road_select-10);
        Relay_y = y_DUE(road_select-10);


        % ��ÿ�������ŵ���CUE�û���ѡ���м̽ڵ�ĸ���
        [SINR_CueToD, SINR_RelayToCue] = judgeSINR_RelaytoD(Relay_x,Relay_y,x_D,y_D,x_CUE,y_CUE,array_fastFading,array_slowFading,bandwidth,Power_UE);
        SINR_CueToD = 10*log10(SINR_CueToD);
        SINR_RelayToCue = 10*log10(SINR_RelayToCue); 

        % �ҵ���������ȵ��ŵ�
        canUsedChannel = intersect(find(SINR_CueToD > SINR_require),find(SINR_RelayToCue> SINR_require));
        max_SINR= max(SINR_CueToD(canUsedChannel));

        if(~isempty(max_SINR))  % ���������������ŵ���ȡ��������
            % ���㴫������
            R = bandwidth * log2(1 + max_SINR); % bit/s
            % ����һ��ʱ϶�ܴ����������
            tranSize = floor(R * (1/1000000));

            if tranSize >= data_save(road_select-10)  % ���Դ�������������ڴ洢��������
                % ���㴫����Ҫ���ѵ�ʱ��
                t_tran = data_save(road_select-10) / R;
                %��������
                total_receive = total_receive + data_save(road_select-10); % �����ܽ��յ�������
                data_save(road_select - 10) = 0;    % ���½ڵ��д洢��������
                road_array(road_select) = 0;
                t_residue(road_select - 10) =  t_residue(road_select - 10) - t_tran*1000000; % ���½ڵ�ʣ�๤��ʱ��           
            else
                 total_receive = total_receive + tranSize; % �����ܽ��յ�������
                 data_save(road_select - 10) = data_save(road_select - 10) - tranSize;    % ���½ڵ��д洢��������
                 t_residue(road_select - 10) =  t_residue(road_select - 10) - 1; % ���½ڵ�ʣ�๤��ʱ��       
            end
        else
            m1 = m1 + 1;
        end
    else  % S -> Relay
        Relay_x = x_DUE(road_select);
        Relay_y = y_DUE(road_select); 

        % ��ÿ�������ŵ���CUE�û���ѡ���м̽ڵ�ĸ���
        [SINR_CueToRelay, SINR_SToCue] = judgeSINR_StoRelay(Relay_x,Relay_y,x_S,y_S,x_CUE,y_CUE,array_fastFading,array_slowFading,bandwidth,Power_UE);
        SINR_CueToRelayDB = 10*log10(SINR_CueToRelay);
        SINR_SToCueDB = 10*log10(SINR_SToCue);

        % �ҵ���������ȵ��ŵ� 
        canUsedChannel = intersect(find(SINR_CueToRelayDB > SINR_require),find(SINR_SToCueDB> SINR_require));
        max_SINR= max(SINR_CueToRelay(canUsedChannel));
        if(~isempty(max_SINR))  % ���������������ŵ���ȡ��������
            % ���㴫������
            R = bandwidth * log2(1 + max_SINR); % bit/s
            % ����һ��ʱ϶�ܴ����������
            tranSize = floor(R * (1/1000000));

            % �ж��м̽ڵ������Ϳռ��Ƿ��㹻
            % ����ڵ�ʣ���������Թ�����ʱ��t_residue s
            % ����ڵ����������ʴ�����Щ��������Ҫ��ʱ��t_tran s
            t_tran = tranSize/R_min;
            if data_save(road_select) + tranSize <= Relay_buffer * M_packet_length && t_residue(road_select) >= t_tran*1000000
                % ����Ӧ�м����������������Ҹ����м̽ڵ�洢��������
                data_save(road_select) = data_save(road_select) + tranSize;
                road_array(road_select + 10) = road_select + 10;
            else
                m2 = m2 + 1;
            end
        else
            m1 = m1 + 1;
        end
    end
   %data_save_temp = [data_save_temp;data_save];
   if total_receive / 102400 > flag
       total_receive
       i
       t_residue
       flag = flag + 1;
   end
    if (~isempty(find(t_residue < 10))) % ��һ���ڵ�����С��10
        t_residue
        total_receive
        i
        m1
        m2
        break;
    end
end
    tranTime = i;
    totalData = total_receive;
    %data_save_time = data_save_temp([2:end],:);

end