clc,clear;
%% ����ͼ����
% ��վ���Ƿ�ΧԲ
alpha = 0:pi/50:2*pi;
R_BS = 500;
x_BS = 0;
y_BS = 0;

x_cir_BS = R_BS*cos(alpha) + x_BS;
y_cir_BS = R_BS*sin(alpha) + x_BS;
plot(x_BS,y_BS,'k^','MarkerFaceColor','k');
hold on;
plot(x_cir_BS,y_cir_BS,'c--');
axis equal;
hold on;

% D2D ͨ�����긲�Ƿ�Χ�����跢�䷽��0,300�������շ���0,480��������180m
x_S = 0;
y_S = 300;
x_D = 0;
y_D = 480;


R_D2D = sqrt((x_D - x_S)^2 +(y_D - y_S)^2 );
x_cir_D2D =  R_D2D*cos(alpha)  + x_S;
y_cir_D2D =  R_D2D*sin(alpha)  + y_S;
plot(x_cir_D2D,y_cir_D2D,'r--');
axis equal;
hold on;

plot(x_S,y_S,'ro','MarkerFaceColor','r');
hold on;
plot(x_D,y_D,'rs','MarkerFaceColor','r');
hold on;
% �������CUE �û� 20��
x = 1000*rand(1,10000) - 500;
y = 1000*rand(1,10000) - 500;
x_CUE_temp = x(x.^2+y.^2<500*500);
y_CUE_temp = y(x.^2+y.^2<500*500);
x_CUE = x_CUE_temp(1:20);
y_CUE = y_CUE_temp(1:20);
scatter(x_CUE,y_CUE,'bx');
hold on;

% �������DUE �û� 10��
x = 2*R_D2D*rand(1,10000) - R_D2D;
y = 2*R_D2D*rand(1,10000) - R_D2D;
x_DUE_temp = x(x.^2+y.^2<R_D2D*R_D2D);
y_DUE_temp = y(x.^2+y.^2<R_D2D*R_D2D);
x_DUE = x_DUE_temp(1:10) + x_S;
y_DUE = y_DUE_temp(1:10) + y_S;
scatter(x_DUE,y_DUE,'mh');
hold on;
axis([-500 500 -500 500]);

%% ��������
M_packet = 50;  %������ 50
M_packet_length = 1024; % ����С 1024bit
Relay_buffer = 10; % �ܴ洢������10
Noise_density = -174; % �����������ܶ� -174dbm/hz    dbm = 10*lg(mw)
Path_loss_exponent = 4; % ·��ϵ�� ��
K = 0.01; % ·���� K
Power_UE = 24; % �û��ķ��书�� 24dbm    10^(24/10) mw;
Energy_loss_factor = 0.3; %��ص�������ʧ30%
T_slot = 1; % һ��ʱ϶ʱ�� 1us   1��=1000����(ms)1��=1000000 ΢��(��s)
SINR_require = 10; % ��Ҫ���������� 10db

Num_CUE = 20; % CUE�û�20��
Num_DUE = 10; % DUE�û�10��
 

% �ɱ����
bandwidth = 720000; % ����720 000 hz = 720khz
R_min = bandwidth*log2(1 + 10);
% �����м̽ڵ������ ������Ϊ��������


%% ����һ ���ѡ����· 
% array_fastFading =  exprnd(1);  % ��ֵΪ1��ָ���ֲ�
% array_slowFading = lognrnd(0,8); %��ֵΪ0������Ϊ8db����̬�����ֲ�


road_array = [1 2 3 4 5 6 7 8 9 10 0 0 0 0 0 0 0 0 0 0]; % ��¼��Щ��·��ѡ��
data_save = [0 0 0 0 0 0 0 0 0 0]; % ��¼�ڵ��еĴ洢�������ĸ���
road_array_temp = find(road_array > 0);
index = randperm(length(road_array_temp));
road_select = road_array_temp(index(1)); % ѡ�����·��

% �������·�������
% 1. ѡ���������ߵ��Ǹ������ŵ�
% 2. ���㴫�����ʣ��Ƿ������������ȣ����������һ������������һ��ʱ϶����ѡ��
% 3. ���㴫���������
% 4. ����ѡ�����ǰ����·���Ǻ�����·���������
% ���1�����ѡ��·����ǰ�����·������ڵ����������������Ƿ��㹻���䣿����ʣ��ռ��Ƿ��㹻��
% �㹻���ʹ�����������һ��ʱ϶����ѡ����·
% ���ѡ����Ǻ�����·������ѡ��ڵ��е����ݡ�

array_fastFading =  exprnd(1,1,100000); % ����ָ���ֲ�����,�������ݣ�ȡ20��
array_slowFading =  lognrnd(0,8,1,100000); % ������̬�����ֲ����飬�������ݣ�ȡ20��

% ��ȡѡ��·����DUE������
if road_select > 10  % Relay->D
    Relay_x = x_DUE(road_select-10);
    Relay_y = y_DUE(road_select-10);
 
    
else  % S -> Relay
    Relay_x = x_DUE(road_select);
    Relay_y = y_DUE(road_select); 
    
    % ��ÿ�������ŵ���CUE�û���ѡ���м̽ڵ�ĸ���
    [SINR_CueToRelay, SINR_RelayToCue] = judgeSINR_StoRelay(Relay_x,Relay_y,x_S,y_S,x_CUE,y_CUE,array_fastFading,array_slowFading,bandwidth,Power_UE);
    SINR_CueToRelayDB = 10*log10(SINR_CueToRelay);
    SINR_RelayToCueDB = 10*log10(SINR_RelayToCue);

    % �ҵ���������ȵ��ŵ� 
    canUsedChannel = intersect(find(SINR_CueToRelayDB > SINR_require),find(SINR_RelayToCueDB> SINR_require));
    max_SINR= max(SINR_CueToRelay(canUsedChannel));
    if(length(max_SINR))  % ���������������ŵ���ȡ��������
        % ���㴫������
        R = bandwidth * log2(1 + max_SINR); % bit/s
        % ����һ��ʱ϶�ܴ����������
        tranSize = R * (1/1000000);
        
        
        
        % �ж��м̽ڵ������Ϳռ��Ƿ��㹻
        % ����ڵ�ʣ���������Թ�����ʱ��t_residue
        % ����ڵ����������ʴ�����Щ��������Ҫ��ʱ��t_tran s
        t_residue = 1111;
        t_tran = tranSize/R_min;
        if data_save(road_select) + tranSize <= Relay_buffer * M_packet && t_residue >= t_tran
            % ����Ӧ�м����������������Ҹ����м̽ڵ�洢��������
            data_save(road_select) = data_save(road_select) + tranSize;
            road_array(road_select + 10) = road_select + 10;
        end
    end
end



 %% ������ ������Ϊ׼��ѡ����·


%% ������ ͬʱ�������ʣ�����������ռ�ѡ����·



