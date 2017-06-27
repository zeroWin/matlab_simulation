% 函数功能 ： 计算节点剩余工作时间
% 输入：
% 1.BatteryCap 节点的电池容量 mAH
% 2.Power_UE  中继的发射功率 dbm
% 3.Energy_loss_factor 能量转换效率
% 输出：last_time 剩余工作时间  s
function last_time = judgeRelayWorkTime(BatteryCap,Power_UE,Energy_loss_factor)
    I = 10^(Power_UE/10) * (1 + Energy_loss_factor) / 


end