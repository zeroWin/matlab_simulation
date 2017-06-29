% �������� �� �����м̽ڵ㣬����R->D���и���CUEƵ�ʵ�SINR������R������CUE��SINR
% ���룺
% 1.ѡ���м̽ڵ������                       Relay_x,Relay_y
% 2.Ŀ�Ľڵ������                             D_x,D_y
% 3.ÿ��CUE�û����� 20��                 CUE_x,CUE_y
% 4.���ɵĿ�˥�����                          array_fastFading
% 5.���ɵ���˥�����                          array_slowFading
% 6.�ŵ�����                                       BandWidth    hz
% 7.�ڵ㷢�书��                                Power_UE      dbm

% ����
% 1.����CUE�û���D��SINR
% 2.Relay������CUE�û���SINR
function [SINR_CueToD, SINR_RelayToCue] = judgeSINR_RelaytoD(Relay_x,Relay_y,D_x,D_y,CUE_x,CUE_y,array_fastFading,array_slowFading,BandWidth,Power_UE)
    K = 0.01; % ·���� K
    a = 4; % ·��ϵ�� ��
    noise_density = -174; % �����ܶ� -174dbm/hz
    P_noise = 10^(noise_density/10) * BandWidth; % ��������
    
    %% ��������CUE�û����м̽ڵ��SINR
    % 1.����Relay->D���书��
    L_RelaytoD = sqrt( (D_x -  Relay_x)^2 + (D_y - Relay_y)^2);
    h_RelaytoD = K * array_fastFading(1) * array_slowFading(1) *  L_RelaytoD^(-a);% �ŵ�ϵ��
    P_RelaytoD = 10^(Power_UE/10) * h_RelaytoD;
    
    % 2.����Cue->D�ĸ��Ź���
    L_CueToD = sqrt( (CUE_x -  D_x).^2 + (CUE_y - D_y).^2);
    h_CueToD = K .* array_fastFading(2:21) .* array_slowFading(2:21) .*  L_CueToD.^(-a);% �ŵ�ϵ��
    P_CueToD = 10^(Power_UE/10).*h_CueToD;
    
    % 3.����SINR
    SINR_CueToD = P_RelaytoD./(P_CueToD + P_noise);
    
    %% ����S������CUE�û���SINR �����ŵ�
    % 1.����CUE->BS���书��
    L_CutToBS = sqrt( (CUE_x -  0).^2 + (CUE_y - 0).^2);
    h_CutToBS = K .* array_fastFading(22:41) .* array_slowFading(22:41) .*  L_CutToBS.^(-a);% �ŵ�ϵ��
    P_CutToBS = 10^(Power_UE/10).*h_CutToBS;
    
    % 2.����Relay->BS����
    L_RelaytoBS = sqrt( (Relay_x - 0)^2 + (Relay_y - 0)^2);
    h_RelaytoBS = K * array_fastFading(42) * array_slowFading(42) *  L_RelaytoBS^(-a);% �ŵ�ϵ��
    P_RelaytoBS = 10^(Power_UE/10) * h_RelaytoBS;
    
    % 3.����SINR
    SINR_RelayToCue = P_CutToBS / (P_RelaytoBS + P_noise);

end