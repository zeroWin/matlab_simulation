% �������� �� ����ڵ�ʣ�๤��ʱ��
% ���룺
% 1.BatteryCap �ڵ�ĵ������ mAH
% 2.Power_UE  �м̵ķ��书�� dbm
% 3.Energy_loss_factor ����ת��Ч��
% 4.�м̽ڵ㹤����ѹ V

% �����last_time ʣ�๤��ʱ��  s
function last_time = judgeRelayWorkTime(BatteryCap,Power_UE,Energy_loss_factor,OperaVol)
    a = 1.3;
    I = 10^(Power_UE/10) * (1 + Energy_loss_factor) / OperaVol; %mA
    last_time = 3600 * BatteryCap /  I^a; % s
end