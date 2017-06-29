% 函数功能 ： 计算发送节点对所有DUE用户的干扰
% 输入：
% 1.发送节点的坐标                             S_x,S_y
% 2.接收节点坐标                                D_x,D_y
% 2.每个CUE用户坐标 20个                 CUE_x,CUE_y
% 3.生成的快衰落矩阵                          array_fastFading      需要长度21
% 4.生成的慢衰落矩阵                          array_slowFading     需要长度21
% 5.信道带宽                                       BandWidth    hz
% 6.节点发射功率                                Power_UE      dbm

% 返回
% 1.所有CUE用户到选择的DUE用户的SINR
function [SINR_CueToD] = judgeSINR_CuetoD( S_x,S_y,D_x,D_y,CUE_x,CUE_y,array_fastFading,array_slowFading,BandWidth,Power_UE)
    K = 0.01; % 路损常数 K
    a = 4; % 路损系数 α
    noise_density = -174; % 噪声密度 -174dbm/hz
    P_noise = 10^(noise_density/10) * BandWidth; % 噪声功率
    
    %% 所有CUE用户到选择的DUE用户的SINR 上行信道
    % 1.计算DUE 发送到接收方的功率
    L_SToD = sqrt( (D_x -  S_x)^2 + (D_y - S_y)^2);
    h_SToD = K * array_fastFading(1) * array_slowFading(1) *  L_SToD^(-a);% 信道系数
    P_SToD = 10^(Power_UE/10)*h_SToD;
    
    % 2.计算CUE->D干扰
    L_CuetoD = sqrt( (D_x - CUE_x).^2 + (D_y - CUE_y).^2);
    h_CuetoD = K * array_fastFading(2:21) .* array_slowFading(2:21) .*  L_CuetoD.^(-a);% 信道系数
    P_CuetoD = 10^(Power_UE/10) * h_CuetoD;
    
    % 3.计算SINR
    SINR_CueToD = P_SToD ./ (P_CuetoD + P_noise);

end