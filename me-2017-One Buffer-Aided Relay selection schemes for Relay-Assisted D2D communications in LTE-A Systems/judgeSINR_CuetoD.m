% �������� �� ���㷢�ͽڵ������DUE�û��ĸ���
% ���룺
% 1.���ͽڵ������                             S_x,S_y
% 2.���սڵ�����                                D_x,D_y
% 2.ÿ��CUE�û����� 20��                 CUE_x,CUE_y
% 3.���ɵĿ�˥�����                          array_fastFading      ��Ҫ����21
% 4.���ɵ���˥�����                          array_slowFading     ��Ҫ����21
% 5.�ŵ�����                                       BandWidth    hz
% 6.�ڵ㷢�书��                                Power_UE      dbm

% ����
% 1.����CUE�û���ѡ���DUE�û���SINR
function [SINR_CueToD] = judgeSINR_CuetoD( S_x,S_y,D_x,D_y,CUE_x,CUE_y,array_fastFading,array_slowFading,BandWidth,Power_UE)
    K = 0.01; % ·���� K
    a = 4; % ·��ϵ�� ��
    noise_density = -174; % �����ܶ� -174dbm/hz
    P_noise = 10^(noise_density/10) * BandWidth; % ��������
    
    %% ����CUE�û���ѡ���DUE�û���SINR �����ŵ�
    % 1.����DUE ���͵����շ��Ĺ���
    L_SToD = sqrt( (D_x -  S_x)^2 + (D_y - S_y)^2);
    h_SToD = K * array_fastFading(1) * array_slowFading(1) *  L_SToD^(-a);% �ŵ�ϵ��
    P_SToD = 10^(Power_UE/10)*h_SToD;
    
    % 2.����CUE->D����
    L_CuetoD = sqrt( (D_x - CUE_x).^2 + (D_y - CUE_y).^2);
    h_CuetoD = K * array_fastFading(2:21) .* array_slowFading(2:21) .*  L_CuetoD.^(-a);% �ŵ�ϵ��
    P_CuetoD = 10^(Power_UE/10) * h_CuetoD;
    
    % 3.����SINR
    SINR_CueToD = P_SToD ./ (P_CuetoD + P_noise);

end