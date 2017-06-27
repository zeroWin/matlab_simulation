% 函数功能 ： 给定中继节点，计算S-R所有复用CUE频率的SINR，计算R到所有CUE的SINR
% 输入：
% 1.选择中继节点的坐标                       Relay_x,Relay_y
% 2.发送节点的坐标                             S_x,S_y
% 3.每个CUE用户坐标 20个                 CUE_x,CUE_y
% 4.生成的快衰落矩阵                          array_fastFading
% 5.生成的慢衰落矩阵                          array_slowFading
% 6.信道带宽                                       BandWidth    hz
% 7.节点发射功率                                Power_UE      dbm

% 返回
% 1.所有CUE用户到Relay的SINR
% 2.Relay到所有CUE用户的SINR
function [SINR_CueToRelay, SINR_RelayToCue] = judgeSINR_StoRelay(Relay_x,Relay_y,S_x,S_y,CUE_x,CUE_y,array_fastFading,array_slowFading,BandWidth,Power_UE)
    K = 0.01; % 路损常数 K
    a = 4; % 路损系数 α
    noise_density = -174; % 噪声密度 -174dbm/hz
    P_noise = 10^(noise_density/10) * BandWidth; % 噪声功率
    
    %% 计算所有CUE用户到中继节点的SINR
    % 1.计算S->Relay发射功率
    L_StoRelay = sqrt( (S_x -  Relay_x)^2 + (S_y - Relay_y)^2);
    h_StoRealy = K * array_fastFading(1) * array_slowFading(1) *  L_StoRelay^(-a);% 信道系数
    P_StoRelay = 10^(Power_UE/10) * h_StoRealy;
    
    % 2.计算Cue->Relay的干扰功率
    L_CueToRelay = sqrt( (CUE_x -  Relay_x).^2 + (CUE_y - Relay_y).^2);
    h_CueToRelay = K .* array_fastFading(2:21) .* array_slowFading(2:21) .*  L_CueToRelay.^(-a);% 信道系数
    P_CueToRelay = 10^(Power_UE/10).*h_CueToRelay;
    
    % 3.计算SINR
    SINR_CueToRelay = P_StoRelay./(P_CueToRelay + P_noise);
    
    %% 计算S到所有CUE用户的SINR 上行信道
    % 1.计算CUE->BS发射功率
    L_CutToBS = sqrt( (CUE_x -  0).^2 + (CUE_y - 0).^2);
    h_CutToBS = K .* array_fastFading(22:41) .* array_slowFading(22:41) .*  L_CutToBS.^(-a);% 信道系数
    P_CutToBS = 10^(Power_UE/10).*h_CutToBS;
    
    % 2.计算S->BS干扰
    L_StoBS = sqrt( (S_x - 0)^2 + (S_y - 0)^2);
    h_StoBS = K * array_fastFading(42) * array_slowFading(42) *  L_StoBS^(-a);% 信道系数
    P_StoBS = 10^(Power_UE/10) * h_StoBS;
    
    % 3.计算SINR
    SINR_RelayToCue = P_CutToBS / (P_StoBS + P_noise);

end