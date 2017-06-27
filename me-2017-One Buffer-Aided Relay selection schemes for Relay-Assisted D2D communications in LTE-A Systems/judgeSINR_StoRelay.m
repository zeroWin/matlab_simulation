% �������� �� �����м̽ڵ㣬����S-R���и���CUEƵ�ʵ�SINR������R������CUE��SINR
% ���룺
% 1.ѡ���м̽ڵ������                       Relay_x,Relay_y
% 2.���ͽڵ������                             S_x,S_y
% 3.ÿ��CUE�û����� 20��                 CUE_x,CUE_y
% 4.���ɵĿ�˥�����                          array_fastFading
% 5.���ɵ���˥�����                          array_slowFading
% 6.�ŵ�����                                       BandWidth    hz
% 7.�ڵ㷢�书��                                Power_UE      dbm

% ����
% 1.����CUE�û���Relay��SINR
% 2.Relay������CUE�û���SINR
function [SINR_CueToRelay, SINR_RelayToCue] = judgeSINR_StoRelay(Relay_x,Relay_y,S_x,S_y,CUE_x,CUE_y,array_fastFading,array_slowFading,BandWidth,Power_UE)
    K = 0.01; % ·���� K
    a = 4; % ·��ϵ�� ��
    noise_density = -174; % �����ܶ� -174dbm/hz
    P_noise = 10^(noise_density/10) * BandWidth; % ��������
    
    %% ��������CUE�û����м̽ڵ��SINR
    % 1.����S->Relay���书��
    L_StoRelay = sqrt( (S_x -  Relay_x)^2 + (S_y - Relay_y)^2);
    h_StoRealy = K * array_fastFading(1) * array_slowFading(1) *  L_StoRelay^(-a);% �ŵ�ϵ��
    P_StoRelay = 10^(Power_UE/10) * h_StoRealy;
    
    % 2.����Cue->Relay�ĸ��Ź���
    L_CueToRelay = sqrt( (CUE_x -  Relay_x).^2 + (CUE_y - Relay_y).^2);
    h_CueToRelay = K .* array_fastFading(2:21) .* array_slowFading(2:21) .*  L_CueToRelay.^(-a);% �ŵ�ϵ��
    P_CueToRelay = 10^(Power_UE/10).*h_CueToRelay;
    
    % 3.����SINR
    SINR_CueToRelay = P_StoRelay./(P_CueToRelay + P_noise);
    
    %% ����S������CUE�û���SINR �����ŵ�
    % 1.����CUE->BS���书��
    L_CutToBS = sqrt( (CUE_x -  0).^2 + (CUE_y - 0).^2);
    h_CutToBS = K .* array_fastFading(22:41) .* array_slowFading(22:41) .*  L_CutToBS.^(-a);% �ŵ�ϵ��
    P_CutToBS = 10^(Power_UE/10).*h_CutToBS;
    
    % 2.����S->BS����
    L_StoBS = sqrt( (S_x - 0)^2 + (S_y - 0)^2);
    h_StoBS = K * array_fastFading(42) * array_slowFading(42) *  L_StoBS^(-a);% �ŵ�ϵ��
    P_StoBS = 10^(Power_UE/10) * h_StoBS;
    
    % 3.����SINR
    SINR_RelayToCue = P_CutToBS / (P_StoBS + P_noise);

end