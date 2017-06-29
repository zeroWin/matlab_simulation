% 函数功能 ： 计算发送节点对所有DUE用户的干扰
% 输入：
% 1.发送节点的坐标                             S_x,S_y
% 2.每个CUE用户坐标 20个                 CUE_x,CUE_y
% 3.生成的快衰落矩阵                          array_fastFading      需要长度21
% 4.生成的慢衰落矩阵                          array_slowFading     需要长度21
% 5.信道带宽                                       BandWidth    hz
% 6.节点发射功率                                Power_UE      dbm

% 返回
% 1.S到所有CUE用户的SINR
function [SINR_SToCue] = judgeSINR_StoCue(S_x,S_y,CUE_x,CUE_y,array_fastFading,array_slowFading,BandWidth,Power_UE)
    K = 0.01; % 路损常数 K
    a = 4; % 路损系数 α
    noise_density = -174; % 噪声密度 -174dbm/hz
    P_noise = 10^(noise_density/10) * BandWidth; % 噪声功率
    
    %% 计算S到所有CUE用户的SINR 上行信道
    % 1.计算CUE->BS发射功率
    L_CutToBS = sqrt( (CUE_x -  0).^2 + (CUE_y - 0).^2);
    h_CutToBS = K .* array_fastFading(1:20) .* array_slowFading(1:20) .*  L_CutToBS.^(-a);% 信道系数
    P_CutToBS = 10^(Power_UE/10).*h_CutToBS;
    
    % 2.计算S->BS干扰
    L_StoBS = sqrt( (S_x - 0)^2 + (S_y - 0)^2);
    h_StoBS = K * array_fastFading(21) * array_slowFading(21) *  L_StoBS^(-a);% 信道系数
    P_StoBS = 10^(Power_UE/10) * h_StoBS;
    
    % 3.计算SINR
    SINR_SToCue = P_CutToBS / (P_StoBS + P_noise);

end