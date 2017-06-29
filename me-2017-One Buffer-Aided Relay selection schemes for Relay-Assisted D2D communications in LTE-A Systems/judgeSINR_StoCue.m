% �������� �� ���㷢�ͽڵ������DUE�û��ĸ���
% ���룺
% 1.���ͽڵ������                             S_x,S_y
% 2.ÿ��CUE�û����� 20��                 CUE_x,CUE_y
% 3.���ɵĿ�˥�����                          array_fastFading      ��Ҫ����21
% 4.���ɵ���˥�����                          array_slowFading     ��Ҫ����21
% 5.�ŵ�����                                       BandWidth    hz
% 6.�ڵ㷢�书��                                Power_UE      dbm

% ����
% 1.S������CUE�û���SINR
function [SINR_SToCue] = judgeSINR_StoCue(S_x,S_y,CUE_x,CUE_y,array_fastFading,array_slowFading,BandWidth,Power_UE)
    K = 0.01; % ·���� K
    a = 4; % ·��ϵ�� ��
    noise_density = -174; % �����ܶ� -174dbm/hz
    P_noise = 10^(noise_density/10) * BandWidth; % ��������
    
    %% ����S������CUE�û���SINR �����ŵ�
    % 1.����CUE->BS���书��
    L_CutToBS = sqrt( (CUE_x -  0).^2 + (CUE_y - 0).^2);
    h_CutToBS = K .* array_fastFading(1:20) .* array_slowFading(1:20) .*  L_CutToBS.^(-a);% �ŵ�ϵ��
    P_CutToBS = 10^(Power_UE/10).*h_CutToBS;
    
    % 2.����S->BS����
    L_StoBS = sqrt( (S_x - 0)^2 + (S_y - 0)^2);
    h_StoBS = K * array_fastFading(21) * array_slowFading(21) *  L_StoBS^(-a);% �ŵ�ϵ��
    P_StoBS = 10^(Power_UE/10) * h_StoBS;
    
    % 3.����SINR
    SINR_SToCue = P_CutToBS / (P_StoBS + P_noise);

end