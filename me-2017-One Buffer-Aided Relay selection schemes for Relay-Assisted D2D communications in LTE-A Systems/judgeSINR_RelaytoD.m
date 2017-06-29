% 函数功能 ： 给定中继节点，计算R->D所有复用CUE频率的SINR，计算R到所有CUE的SINR
% 输入：
% 1.选择中继节点的坐标                       Relay_x,Relay_y
% 2.目的节点的坐标                             D_x,D_y
% 3.每个CUE用户坐标 20个                 CUE_x,CUE_y
% 4.生成的快衰落矩阵                          array_fastFading
% 5.生成的慢衰落矩阵                          array_slowFading
% 6.信道带宽                                       BandWidth    hz
% 7.节点发射功率                                Power_UE      dbm

% 返回
% 1.所有CUE用户到D的SINR
% 2.Relay到所有CUE用户的SINR
function [SINR_CueToD, SINR_RelayToCue] = judgeSINR_RelaytoD(Relay_x,Relay_y,D_x,D_y,CUE_x,CUE_y,array_fastFading,array_slowFading,BandWidth,Power_UE)
    K = 0.01; % 路损常数 K
    a = 4; % 路损系数 α
    noise_density = -174; % 噪声密度 -174dbm/hz
    P_noise = 10^(noise_density/10) * BandWidth; % 噪声功率
    
    %% 计算所有CUE用户到中继节点的SINR
    % 1.计算Relay->D发射功率
    L_RelaytoD = sqrt( (D_x -  Relay_x)^2 + (D_y - Relay_y)^2);
    h_RelaytoD = K * array_fastFading(1) * array_slowFading(1) *  L_RelaytoD^(-a);% 信道系数
    P_RelaytoD = 10^(Power_UE/10) * h_RelaytoD;
    
    % 2.计算Cue->D的干扰功率
    L_CueToD = sqrt( (CUE_x -  D_x).^2 + (CUE_y - D_y).^2);
    h_CueToD = K .* array_fastFading(2:21) .* array_slowFading(2:21) .*  L_CueToD.^(-a);% 信道系数
    P_CueToD = 10^(Power_UE/10).*h_CueToD;
    
    % 3.计算SINR
    SINR_CueToD = P_RelaytoD./(P_CueToD + P_noise);
    
    %% 计算S到所有CUE用户的SINR 上行信道
    % 1.计算CUE->BS发射功率
    L_CutToBS = sqrt( (CUE_x -  0).^2 + (CUE_y - 0).^2);
    h_CutToBS = K .* array_fastFading(22:41) .* array_slowFading(22:41) .*  L_CutToBS.^(-a);% 信道系数
    P_CutToBS = 10^(Power_UE/10).*h_CutToBS;
    
    % 2.计算Relay->BS干扰
    L_RelaytoBS = sqrt( (Relay_x - 0)^2 + (Relay_y - 0)^2);
    h_RelaytoBS = K * array_fastFading(42) * array_slowFading(42) *  L_RelaytoBS^(-a);% 信道系数
    P_RelaytoBS = 10^(Power_UE/10) * h_RelaytoBS;
    
    % 3.计算SINR
    SINR_RelayToCue = P_CutToBS / (P_RelaytoBS + P_noise);

end