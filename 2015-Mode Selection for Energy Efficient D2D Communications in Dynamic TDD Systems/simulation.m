clc;clear;
T = 1;
p10_max = 0.25;
p12_max = 0.25;
p02_max = 40;
G0 = 5.7 * 10^(-4);
a = 4;
W = 5 * 10^6; %hz
CF = 1 * 10^9; %hz
noise = (10^(-17.4)) * W / 1000;
Rcell = 500;

%% D10 = 250m k = 1 mobile user energy consumption ok
D10 = 250;
k = 1;
r10_max = W*log(1 + (p10_max * G0 * D10^(-a))/ noise);
r02_max = W*log(1 + (p02_max * G0 * Rcell^(-a))/ noise);
btgt = k * T * r10_max * r02_max / (r10_max + r02_max);
b = (exp(btgt/(W*T)) - 1) * noise;
resultD12 = (p12_max * G0 / b)^(1/a); % User-1 power feasibility radius

% »­Í¼
h = figure(1);
set(h,'name','D2D optimality area when mobile user energy consumption','Numbertitle','off')
subplot(1,2,1);
x_user1 = 0;
y_user1 = 250;

%plot bs user-1 power feasible and user-1 energy-optimal
t = linspace(0,2*pi);
xBS = Rcell * cos(t);
yBS = Rcell * sin(t);
plot(xBS,yBS,'k');
hold on;
xUser1 = x_user1 + resultD12*cos(t);
yUser1 = y_user1 + resultD12*sin(t);
f = fill(xUser1,yUser1,[111/255,204/255,221/255]);
set(f,'EdgeColor','k','FaceAlpha',0.5);
axis equal; 
% D2D more energy efficient
r_all = linspace(50,500,30);
for i=1:1:30
    r_user2 = r_all(i); %User-2 location
    t = linspace(0,2*pi,500);
    x_user2 = r_user2*cos(t);
    y_user2 = r_user2*sin(t);
    D12 = (x_user2.^2 + (y_user2 - D10).^2).^(1/2);
    D02 = r_user2;
    r02_max_temp = W*log(1+p02_max*G0*(D02^(-a))/noise);
    t10_xin = T - btgt/r02_max_temp;
    E_D2D = (exp(btgt/(W*T))-1)*T./(G0*D12.^(-a));
    E_10 = (exp(btgt/(W*t10_xin))-1)*t10_xin/(G0*D10^(-a));
    
    get_location = E_D2D < E_10;
    plot(x_user2(get_location),y_user2(get_location),'r.');
    hold on;
end
text(10,-30,'BS','fontsize',10)
text(0,y_user1-30,'User-1','fontsize',10)
hold on
plot (0,0,'.k')
plot (x_user1,y_user1,'.k')
title('k = 1');
xlabel('(a) D2D optimality area for D10 = 250m');  %xÖá


%% D10 = 450m k = 1 mobile user energy consumption ok
D10 = 450;
k = 1;
r10_max = W*log(1 + (p10_max * G0 * D10^(-a))/ noise);
r02_max = W*log(1 + (p02_max * G0 * Rcell^(-a))/ noise);
btgt = k * T * r10_max * r02_max / (r10_max + r02_max);
b = (exp(btgt/(W*T)) - 1) * noise;
resultD12 = (p12_max * G0 / b)^(1/a);

% »­Í¼
x_user1 = 0;
y_user1 = 450;
subplot(1,2,2);
%plot bs user-1 power feasible and user-1 energy-optimal
t = linspace(0,2*pi);
xBS = Rcell * cos(t);
yBS = Rcell * sin(t);
plot(xBS,yBS,'k');
hold on;
xUser1 = x_user1 + resultD12*cos(t);
yUser1 = y_user1 + resultD12*sin(t);
f = fill(xUser1,yUser1,[111/255,204/255,221/255]);
set(f,'EdgeColor','k','FaceAlpha',0.5);
axis equal;
% D2D more energy efficient
r_all = linspace(50,500,30);
for i=1:1:30
    r_user2 = r_all(i); %User-2 location
    t = linspace(0,2*pi,500);
    x_user2 = r_user2*cos(t);
    y_user2 = r_user2*sin(t);
    D12 = (x_user2.^2 + (y_user2 - D10).^2).^(1/2);
    D02 = r_user2;
    r02_max_temp = W*log(1+p02_max*G0*(D02^(-a))/noise);
    t10_xin = T - btgt/r02_max_temp;
    E_D2D = (exp(btgt/(W*T))-1)*T./(G0*D12.^(-a));
    E_10 = (exp(btgt/(W*t10_xin))-1)*t10_xin/(G0*D10^(-a));
    
    get_location = E_D2D < E_10;
    plot(x_user2(get_location),y_user2(get_location),'r.');
    hold on;
end
text(10,-30,'BS','fontsize',10)
text(0,y_user1-30,'User-1','fontsize',10)
hold on
plot (0,0,'.k')
plot (x_user1,y_user1,'.k')
title('k = 1');
xlabel('(b) D2D optimality area for D10 = 450m');  %xÖá

%% D10 = 250m k = 0.5 overall network energy consumption
D10 = 250;
k = 0.5;
r10_max = W*log(1 + (p10_max * G0 * D10^(-a))/ noise);
r02_max = W*log(1 + (p02_max * G0 * Rcell^(-a))/ noise);

btgt = k * T * r10_max * r02_max / (r10_max + r02_max);
b = (exp(btgt/(W*T)) - 1) * noise;
resultD12 = (p12_max * G0 / b)^(1/a);


% »­Í¼
x_user1 = 0;
y_user1 = 250;
h = figure(2);
set(h,'name','D2D optimality area when minimizing the network energy consumption,D10 = 250m, D02 ranges from 50m to 500m, and k = 0.5','Numbertitle','off')
subplot(1,2,1);

%plot bs user-1 power feasible and user-1 energy-optimal
t = linspace(0,2*pi);
xBS = Rcell * cos(t);
yBS = Rcell * sin(t);
plot(xBS,yBS,'k');
hold on;
xUser1 = x_user1 + resultD12*cos(t);
yUser1 = y_user1 + resultD12*sin(t);
f = fill(xUser1,yUser1,[111/255,204/255,221/255]);
set(f,'EdgeColor','k','FaceAlpha',0.5);
axis equal;

% D2D more energy efficient
r_all = linspace(50,500,30);
for i=1:1:30
    r_user2 = r_all(i); %User-2 location
    t = linspace(0,2*pi,500);
    x_user2 = r_user2*cos(t);
    y_user2 = r_user2*sin(t);
    D12 = (x_user2.^2 + (y_user2 - D10).^2).^(1/2);
    D02 = (x_user2(1)^2 + y_user2(1)^2)^(1/2);
    r02_max_temp = W*log(1+p02_max*G0*(D02^(-a))/noise);
    r10_max_temp = r10_max;
    
    t02_min = btgt/r02_max_temp;
    t02_max = T-btgt/r10_max_temp;
    % D2D power
    E_D2D = (exp(btgt/(W*T))-1)*T./(G0*D12.^(-a));
    
    % All network power
    G10 = G0*D10^(-a);
    G02 = G0*D02^(-a);
    
    options = optimset('Algorithm','sqp','display','none');
    [x,E_10And02] = fmincon(@(x) findMinAllNetworkEnergy(x,btgt,T,W,G10,G02),t02_min,[],[],[],[],t02_min,t02_max,[],options);
    get_location = E_D2D < E_10And02;
    plot(x_user2(get_location),y_user2(get_location),'r.');
    hold on;
end

text(10,-30,'BS','fontsize',10)
text(0,y_user1-30,'User-1','fontsize',10)
hold on
plot (0,0,'.k')
plot (x_user1,y_user1,'.k')
xlabel('(a) D2D optimality area and D2D power feasibility area');  %xÖá


subplot(1,2,2);
% D2D more energy efficient
r_all = linspace(50,500,30);
for i=1:1:30
    r_user2 = r_all(i); %User-2 location
    t = linspace(0,2*pi,500);
    x_user2 = r_user2*cos(t);
    y_user2 = r_user2*sin(t);
    D12 = (x_user2.^2 + (y_user2 - D10).^2).^(1/2);
    D02 = (x_user2(1)^2 + y_user2(1)^2)^(1/2);
    r02_max_temp = W*log(1+p02_max*G0*(D02^(-a))/noise);
    r10_max_temp = r10_max;
    
    t02_min = btgt/r02_max_temp;
    t02_max = T-btgt/r10_max_temp;
    % D2D power
    E_D2D = (exp(btgt/(W*T))-1)*T./(G0*D12.^(-a));
    
    % All network power
    G10 = G0*D10^(-a);
    G02 = G0*D02^(-a);
    
    options = optimset('Algorithm','sqp','display','none');
    [x,E_10And02] = fmincon(@(x) findMinAllNetworkEnergy(x,btgt,T,W,G10,G02),t02_min,[],[],[],[],t02_min,t02_max,[],options);
    PowerSavePercent = ((E_10And02 - E_D2D)/E_10And02)*100;
    axis([-500 500 -500 500 -1000 500]); 
    plot3(x_user2,y_user2,PowerSavePercent);
    hold on;
end
plot3(0,0,0,'r.','MarkerSize',20);
plot3(0,D10,0,'k.','MarkerSize',20);
text(0,-30,0,'BS','fontsize',10)
text(0,D10-30,0,'User-1','fontsize',10)
title('Energy saving with D2D-mode(%)')
xlabel('(b) Percentage of energt saved in D2D mode compared with cellular mode');  %xÖá

%% D10 = 450m k = 0.5 overall network energy consumption
D10 = 450;
k = 0.5;
r10_max = W*log(1 + (p10_max * G0 * D10^(-a))/ noise);
r02_max = W*log(1 + (p02_max * G0 * Rcell^(-a))/ noise);

btgt = k * T * r10_max * r02_max / (r10_max + r02_max);
b = (exp(btgt/(W*T)) - 1) * noise;
resultD12 = (p12_max * G0 / b)^(1/a);

% »­Í¼
x_user1 = 0;
y_user1 = 450;
h = figure(3);
set(h,'name','D2D optimality area when minimizing the network energy consumption,D10 = 450m, D02 ranges from 50m to 500m, and k = 0.5','Numbertitle','off')
subplot(1,2,1);

%plot bs user-1 power feasible and user-1 energy-optimal
t = linspace(0,2*pi);
xBS = Rcell * cos(t);
yBS = Rcell * sin(t);
plot(xBS,yBS,'k');
hold on;
xUser1 = x_user1 + resultD12*cos(t);
yUser1 = y_user1 + resultD12*sin(t);
f = fill(xUser1,yUser1,[111/255,204/255,221/255]);
set(f,'EdgeColor','k','FaceAlpha',0.5);
axis equal;
% D2D more energy efficient
r_all = linspace(50,500,30);
for i=1:1:30
    r_user2 = r_all(i); %User-2 location
    t = linspace(0,2*pi,500);
    x_user2 = r_user2*cos(t);
    y_user2 = r_user2*sin(t);
    D12 = (x_user2.^2 + (y_user2 - D10).^2).^(1/2);
    D02 = (x_user2(1)^2 + y_user2(1)^2)^(1/2);
    r02_max_temp = W*log(1+p02_max*G0*(D02^(-a))/noise);
    r10_max_temp = r10_max;
    
    t02_min = btgt/r02_max_temp;
    t02_max = T-btgt/r10_max_temp;
    % D2D power
    E_D2D = (exp(btgt/(W*T))-1)*T./(G0*D12.^(-a));
    
    % All network power
    G10 = G0*D10^(-a);
    G02 = G0*D02^(-a);
    
    options = optimset('Algorithm','sqp','display','none');
    [x,E_10And02] = fmincon(@(x) findMinAllNetworkEnergy(x,btgt,T,W,G10,G02),t02_min,[],[],[],[],t02_min,t02_max,[],options);
    get_location = E_D2D < E_10And02;
    plot(x_user2(get_location),y_user2(get_location),'r.');
    hold on;
end

text(10,-30,'BS','fontsize',10)
text(0,y_user1-30,'User-1','fontsize',10)
hold on
plot (0,0,'.k')
plot (x_user1,y_user1,'.k')
xlabel('(a) D2D optimality area and D2D power feasibility area');  %xÖá

subplot(1,2,2);
% D2D more energy efficient
r_all = linspace(50,500,30);
for i=1:1:30
    r_user2 = r_all(i); %User-2 location
    t = linspace(0,2*pi,500);
    x_user2 = r_user2*cos(t);
    y_user2 = r_user2*sin(t);
    D12 = (x_user2.^2 + (y_user2 - D10).^2).^(1/2);
    D02 = (x_user2(1)^2 + y_user2(1)^2)^(1/2);
    r02_max_temp = W*log(1+p02_max*G0*(D02^(-a))/noise);
    r10_max_temp = r10_max;
    
    t02_min = btgt/r02_max_temp;
    t02_max = T-btgt/r10_max_temp;
    % D2D power
    E_D2D = (exp(btgt/(W*T))-1)*T./(G0*D12.^(-a));
    
    % All network power
    G10 = G0*D10^(-a);
    G02 = G0*D02^(-a);
    
    options = optimset('Algorithm','sqp','display','none');
    [x,E_10And02] = fmincon(@(x) findMinAllNetworkEnergy(x,btgt,T,W,G10,G02),t02_min,[],[],[],[],t02_min,t02_max,[],options);
    PowerSavePercent = ((E_10And02 - E_D2D)/E_10And02)*100;
    axis([-500 500 -500 500 -1000 500]); 
    plot3(x_user2,y_user2,PowerSavePercent);
    hold on;
end
plot3(0,0,0,'r.','MarkerSize',20);
plot3(0,D10,0,'k.','MarkerSize',20);
text(0,-30,0,'BS','fontsize',10)
text(0,D10-30,0,'User-1','fontsize',10)
title('Energy saving with D2D-mode(%)')
xlabel('(b) Percentage of energt saved in D2D mode compared with cellular mode');  %xÖá


%% D10 = 250m k = 0.1 0.5 1 overall network energy consumption ok
D10 = 250;


% default
x_user1 = 0;
y_user1 = 250;
h = figure(4);
set(h,'name','D2D optimality area with different target rate,D10 = 250m, D02 ranges from 50m to 500m','Numbertitle','off')
%k = 0.1
k = 0.1;
r10_max = W*log(1 + (p10_max * G0 * D10^(-a))/ noise);
r02_max = W*log(1 + (p02_max * G0 * Rcell^(-a))/ noise);

btgt = k * T * r10_max * r02_max / (r10_max + r02_max);
b = (exp(btgt/(W*T)) - 1) * noise;
resultD12 = (p12_max * G0 / b)^(1/a);
subplot(1,3,1);

%plot bs user-1 power feasible and user-1 energy-optimal
t = linspace(0,2*pi);
xBS = Rcell * cos(t);
yBS = Rcell * sin(t);
plot(xBS,yBS,'k');
hold on;
xUser1 = x_user1 + resultD12*cos(t);
yUser1 = y_user1 + resultD12*sin(t);
f = fill(xUser1,yUser1,[111/255,204/255,221/255]);
set(f,'EdgeColor','k','FaceAlpha',0.5);
axis equal;
% D2D more energy efficient
r_all = linspace(50,500,30);
for i=1:1:30
    r_user2 = r_all(i); %User-2 location
    t = linspace(0,2*pi,500);
    x_user2 = r_user2*cos(t);
    y_user2 = r_user2*sin(t);
    D12 = (x_user2.^2 + (y_user2 - D10).^2).^(1/2);
    D02 = (x_user2(1)^2 + y_user2(1)^2)^(1/2);
    r02_max_temp = W*log(1+p02_max*G0*(D02^(-a))/noise);
    r10_max_temp = r10_max;
    
    t02_min = btgt/r02_max_temp;
    t02_max = T-btgt/r10_max_temp;
    % D2D power
    E_D2D = (exp(btgt/(W*T))-1)*T./(G0*D12.^(-a));
    
    % All network power
    G10 = G0*D10^(-a);
    G02 = G0*D02^(-a);
    
    options = optimset('Algorithm','sqp','display','none');
    [x,E_10And02] = fmincon(@(x) findMinAllNetworkEnergy(x,btgt,T,W,G10,G02),t02_min,[],[],[],[],t02_min,t02_max,[],options);
    get_location = E_D2D < E_10And02;
    plot(x_user2(get_location),y_user2(get_location),'r.');
    hold on;
end

text(10,-30,'BS','fontsize',10)
text(0,y_user1-30,'User-1','fontsize',10)
hold on
plot (0,0,'.k')
plot (x_user1,y_user1,'.k')
xlabel('(a) Rate target is equal to 10% of the maximal feasible rate');  %xÖá

%k = 0.5
k = 0.5;
r10_max = W*log(1 + (p10_max * G0 * D10^(-a))/ noise);
r02_max = W*log(1 + (p02_max * G0 * Rcell^(-a))/ noise);

btgt = k * T * r10_max * r02_max / (r10_max + r02_max);
b = (exp(btgt/(W*T)) - 1) * noise;
resultD12 = (p12_max * G0 / b)^(1/a);
subplot(1,3,2);

%plot bs user-1 power feasible and user-1 energy-optimal
t = linspace(0,2*pi);
xBS = Rcell * cos(t);
yBS = Rcell * sin(t);
plot(xBS,yBS,'k');
hold on;
xUser1 = x_user1 + resultD12*cos(t);
yUser1 = y_user1 + resultD12*sin(t);
f = fill(xUser1,yUser1,[111/255,204/255,221/255]);
set(f,'EdgeColor','k','FaceAlpha',0.5);
axis equal;
% D2D more energy efficient
r_all = linspace(50,500,30);
for i=1:1:30
    r_user2 = r_all(i); %User-2 location
    t = linspace(0,2*pi,500);
    x_user2 = r_user2*cos(t);
    y_user2 = r_user2*sin(t);
    D12 = (x_user2.^2 + (y_user2 - D10).^2).^(1/2);
    D02 = (x_user2(1)^2 + y_user2(1)^2)^(1/2);
    r02_max_temp = W*log(1+p02_max*G0*(D02^(-a))/noise);
    r10_max_temp = r10_max;
    
    t02_min = btgt/r02_max_temp;
    t02_max = T-btgt/r10_max_temp;
    % D2D power
    E_D2D = (exp(btgt/(W*T))-1)*T./(G0*D12.^(-a));
    
    % All network power
    G10 = G0*D10^(-a);
    G02 = G0*D02^(-a);
    
    options = optimset('Algorithm','sqp','display','none');
    [x,E_10And02] = fmincon(@(x) findMinAllNetworkEnergy(x,btgt,T,W,G10,G02),t02_min,[],[],[],[],t02_min,t02_max,[],options);
    get_location = E_D2D < E_10And02;
    plot(x_user2(get_location),y_user2(get_location),'r.');
    hold on;
end

text(10,-30,'BS','fontsize',10)
text(0,y_user1-30,'User-1','fontsize',10)
hold on
plot (0,0,'.k')
plot (x_user1,y_user1,'.k')
xlabel('(b) Rate target is equal to 50% of the maximal feasible rate');  %xÖá

%k = 1
k = 1;
r10_max = W*log(1 + (p10_max * G0 * D10^(-a))/ noise);
r02_max = W*log(1 + (p02_max * G0 * Rcell^(-a))/ noise);

btgt = k * T * r10_max * r02_max / (r10_max + r02_max);
b = (exp(btgt/(W*T)) - 1) * noise;
resultD12 = (p12_max * G0 / b)^(1/a);
subplot(1,3,3);

%plot bs user-1 power feasible and user-1 energy-optimal
t = linspace(0,2*pi);
xBS = Rcell * cos(t);
yBS = Rcell * sin(t);
plot(xBS,yBS,'k');
hold on;
xUser1 = x_user1 + resultD12*cos(t);
yUser1 = y_user1 + resultD12*sin(t);
f = fill(xUser1,yUser1,[111/255,204/255,221/255]);
set(f,'EdgeColor','k','FaceAlpha',0.5);
axis equal;
% D2D more energy efficient
r_all = linspace(50,500,30);
for i=1:1:30
    r_user2 = r_all(i); %User-2 location
    t = linspace(0,2*pi,500);
    x_user2 = r_user2*cos(t);
    y_user2 = r_user2*sin(t);
    D12 = (x_user2.^2 + (y_user2 - D10).^2).^(1/2);
    D02 = (x_user2(1)^2 + y_user2(1)^2)^(1/2);
    r02_max_temp = W*log(1+p02_max*G0*(D02^(-a))/noise);
    r10_max_temp = r10_max;
    
    t02_min = btgt/r02_max_temp;
    t02_max = T-btgt/r10_max_temp;
    % D2D power
    E_D2D = (exp(btgt/(W*T))-1)*T./(G0*D12.^(-a));
    
    % All network power
    G10 = G0*D10^(-a);
    G02 = G0*D02^(-a);
    
    options = optimset('display','none');
    [x,E_10And02] = fmincon(@(x) findMinAllNetworkEnergy(x,btgt,T,W,G10,G02),t02_min,[],[],[],[],t02_min,t02_max,[],options);
    get_location = E_D2D < E_10And02;
    plot(x_user2(get_location),y_user2(get_location),'r.');
    hold on;
end

text(10,-30,'BS','fontsize',10)
text(0,y_user1-30,'User-1','fontsize',10)
hold on
plot (0,0,'.k')
plot (x_user1,y_user1,'.k')
xlabel('(c) Rate target is equal to maximal feasible rate');  %xÖá