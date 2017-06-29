% 函数功能 ： 计算节点剩余工作时间
% 输入：
% 1.BatteryCap 节点的电池容量 mAH
% 2.Power_UE  中继的发射功率 dbm
% 3.Energy_loss_factor 能量转换效率
% 4.中继节点工作电压 V

% 输出：last_time 剩余工作时间  s
function last_time = judgeRelayWorkTime(BatteryCap,Power_UE,Energy_loss_factor,OperaVol)
    a = 1.3;
    I = 10^(Power_UE/10) * (1 + Energy_loss_factor) / OperaVol; %mA
    last_time = 3600 * BatteryCap /  I^a; % s
end