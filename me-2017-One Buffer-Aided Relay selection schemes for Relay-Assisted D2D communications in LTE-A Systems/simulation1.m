clc,clear;
% ��ԭʼ���ĵķ���

%% ��������
Packet_num = 50;  % ���͵İ��ĸ���
Relay_num = 4;    % �м̵ĸ���
P_max = 10;         % mw �����ܹ��ﵽ�������
L_buffer = 10;      % �м��ܹ��洢�ĸ���
Noise = 1;            % mw ��������
SNR = 10;            % Ҫ����������

%% ����1  �̶��м̽ڵ������
% �仯���ͽڵ����������������ͼ
% ͼ1���ܹ����Ͷ��ٸ���
% ͼ2��������50��������Ҫ��ʱ��
E_relay = 600;      % mj �м̽ڵ�����

% ���ܹ����Ͷ��ٸ���

    
for E_source = 100:50:800 % ������100��800
    total_packet = 0;
    total_time = 0;
    for simulationTrials = 1:1:1000 % ����1000��
        Packect_s = Packet_num;
        Packect_d = 0;
        send_error = 0;
        Relay_contain = zeros(1,Relay_num);  % ��ʼ���м̴洢�ж��ٸ�����
        E_relay_con = zeros(1,Relay_num) +   E_relay; %  ��ʼ���м�ʣ������
        E_source_con = E_source;  % ���ͽڵ�ʣ������


        for testTime = 1:1:1000 %����ѡ��1000����·

            % ���ɸ�ʱ����ŵ�ϵ��
            h = normrnd(0,1,1,Relay_num*2);

            % ����s->r ��Ȩ��
            h_sr = h(1:Relay_num).^2; % s->r���ŵ�ϵ��
            P_sm = SNR*Noise./h_sr;     % ÿ�����͵�ÿ���м�����Ҫ����С����
            W_sr = zeros(1,Relay_num);  % ����ÿһ����·��Ȩ��
            e_m =    L_buffer -  Relay_contain; % ����ÿ���м���tʱ�̻��ܴ���ٸ��ڵ�
            s_m =    floor(E_relay_con ./ P_max) -Relay_contain;% ����ÿ���м���tʱ�̻��ܷ��Ͷ��ٸ��ڵ�

            for i = 1:1:Relay_num % ����ÿһ����·��Ȩ��
                if P_sm(i) <= min(P_max, E_source_con)  &&  Packect_s > 0
                    W_sr(i) =  E_source_con/P_sm(i) * min(e_m(i),s_m(i));  
                end
            end


            % ����r->d ��Ȩ��
            h_rd = h(Relay_num+1:end).^2; % r->d���ŵ�ϵ��
            P_rd = SNR*Noise./h_rd;     % ÿ���м̷�������Ҫ����С���� 
            W_rd = zeros(1,Relay_num);  % ����ÿһ����·��Ȩ��

            for i = 1:1:Relay_num % ����ÿһ����·��Ȩ��
                if P_rd(i) <= min(P_max,E_relay_con(i))
                    W_rd(i) = E_relay_con(i) / P_rd(i) * Relay_contain(i);
                end
            end

            [maxNum, maxLocation] = max([W_sr,W_rd]);
            if(maxNum == 0) % ����Ȩ��Ϊ0 ��û��������������· ����50��û��������������·����Ϊ�������㣬�޷���ɷ���
                send_error = send_error + 1;
            else
                if(maxLocation <= Relay_num) % ���ֵ��ǰ�����·
                     Relay_contain(maxLocation) = Relay_contain(maxLocation) + 1;
                     E_source_con  = E_source_con - P_sm(maxLocation); % �޸�sourceʣ������
                     Packect_s = Packect_s - 1;
                else % ���ֵ�ں������·
                     maxLocation = maxLocation - Relay_num;
                     Relay_contain(maxLocation) = Relay_contain(maxLocation) - 1;
                     E_relay_con(maxLocation) = E_relay_con(maxLocation) -  P_rd(maxLocation);
                     Packect_d = Packect_d + 1;
                end        
                send_error = 0;
            end
            if(send_error == 50)
                %s = sprintf('�����޷���ɣ�����Ͱ��ɹ����ո���Ϊ:%d\n,ѭ������Ϊ:%d\n',Packect_d,testTime); 
                %disp(s);
                break;
            end
            if(Packect_d == Packet_num)
               % s = sprintf('���ͳɹ�ѭ������Ϊ:%d\n',testTime); 
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
    s = sprintf('ƽ�����ո���Ϊ:%d\n',floor(total_packet/1000)); 
     disp(s);
    s = sprintf('ƽ�����պķ�ʱ��Ϊ:%d\n',floor(total_time/1000)); 
     disp(s);
end

    




