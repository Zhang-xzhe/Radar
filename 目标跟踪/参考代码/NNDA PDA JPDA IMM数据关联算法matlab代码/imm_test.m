
% 基于IMM算法的目标跟踪
% 
function imm_test()
clear all; 
clc; 
echo off; 
 
%=============================== 建立模型 ===============================
% 仿真参数
simTime=100;      %仿真迭代次数 
T=1;                     %采样时间 
w2=3*2*pi/360;     %模型2转弯率3度 
w3=-3*2*pi/360;    %模型3转弯率-3度 
H=[1,0,0,0;0,0,1,0];                      %模型量测矩阵 
G=[T^2/2,0;T,0;0,T^2/2;0,T];              %模型过程噪声加权矩阵 
r=200;                                 %20 2000 
R=[r,0;0,r];                            %模型量测噪声协方差矩阵 
Q=[10,0;0,10];                                  %模型过程噪声协方差矩阵 
 
F1=[1,T,0,0;0,1,0,0;0,0,1,T;0,0,0,1];     %模型1状态转移矩阵 
 
F2=[1,sin(w2*T)/w2,0,(cos(w2*T)-1)/w2; 
    0,cos(w2*T),0,sin(w2*T); 
    0,(1-cos(w2*T))/w2,1,sin(w2*T)/w2; 
    0,-sin(w2*T),0,cos(w2*T)];            %模型2状态转移矩阵 左转弯 
 
F3=[1,sin(w3*T)/w3,0,(cos(w3*T)-1)/w3; 
    0,cos(w3*T),0,sin(w3*T); 
    0,(1-cos(w3*T))/w3,1,sin(w3*T)/w3; 
    0,-sin(w3*T),0,cos(w3*T)];            %模型3状态转移矩阵 右转弯 
 
w4=3*2*pi/360;     %模型2转弯率3度 
w5=-3*2*pi/360;    %模型3转弯率-3度 
 
F4=[1,sin(w4*T)/w4,0,(cos(w4*T)-1)/w4; 
    0,cos(w4*T),0,sin(w4*T); 
    0,(1-cos(w4*T))/w4,1,sin(w4*T)/w4; 
    0,-sin(w4*T),0,cos(w4*T)];            %模型2状态转移矩阵 
 
F5=[1,sin(w5*T)/w5,0,(cos(w5*T)-1)/w5; 
    0,cos(w5*T),0,sin(w5*T); 
    0,(1-cos(w5*T))/w5,1,sin(w5*T)/w5; 
    0,-sin(w5*T),0,cos(w5*T)];            %模型3状态转移矩阵 
 
x0=[1000,200,1000,200]';  % 初始状态  
% 产生量测数据
%[z,z_true]=targets();
randn('state',sum(100*clock)); % Shuffle the pack! 
x = zeros(4,simTime); 
z = zeros(2,simTime);         %含噪声量测数据 
z_true = zeros(2,simTime);    %真值数据 
measureNoise = zeros(2,simTime); 
measureNoise = sqrt(R)*randn(2,simTime);  %产生量测噪声 
x(:,1)=x0; 
z(:,1)=H*x(:,1)+measureNoise(:,1); 
z_true(:,1)=H*x(:,1); 
for a=2:simTime 
    if (a>=20)&&(a<=40)  
        x(:,a)=F4*x(:,a-1);      %20--40s左转  
    elseif (a>=60)&&(a<=80)  
       x(:,a)=F5*x(:,a-1);        %60--80s右转  
    else 
        x(:,a)=F1*x(:,a-1);      %匀速直线运动 
    end;  
z(:,a)=H*x(:,a)+measureNoise(:,a); 
z_true(:,a)=H*x(:,a); 
end; 
%===============================
%     IMM
%===============================
 
%初始化
x1_IMM = zeros(4,1);      %模型1IMM算法状态估计值 
x2_IMM = zeros(4,1);      %模型2IMM算法状态估计值 
x3_IMM = zeros(4,1);      %模型3IMM算法状态估计值 
x_pro_IMM = zeros(4,simTime);   %IMM算法模型综合状态估计值 
P_IMM=zeros(4,4,simTime);       %IMM算法模型综合状态协方差矩阵 
P1_IMM=zeros(4,4); 
P2_IMM=zeros(4,4); 
P3_IMM=zeros(4,4);              %IMM算法各模型协方差矩阵 
r1_IMM=zeros(2,1); 
r2_IMM=zeros(2,1); 
r3_IMM=zeros(2,1); 
S1_IMM=zeros(2,2); 
S2_IMM=zeros(2,2); 
S3_IMM=zeros(2,2); 
 
 
%初始化
x_pro_IMM(:,1)=x0; 
 
pij=[0.9,0.05,0.05; 
    0.1,0.8,0.1;     0.05,0.15,0.8];    %模型转移概率矩阵 
%pij=[0.6,0.2,0.2;0.2,0.6,0.2;0.25,0.15,0.6];    %模型转移概率矩阵
u_IMM=zeros(3,simTime); 
u_IMM(:,1)=[0.3,0.3,0.4]';  %IMM算法模型概率 
 
x1_IMM=x0;x2_IMM=x0;x3_IMM=x0;  %IMM算法各模型初始状态 
 
P0=diag([1000,500,1000,500]);  %初始状态协方差矩阵 
 
P1_IMM=P0;P2_IMM=P0;P3_IMM=P0; 
 
P_IMM(:,:,1)=P0; 
 
%main loop
 
for t=1:simTime-1 
     
    %第一步Interacting（只针对IMM算法）
    c_j=pij'*u_IMM(:,t); 
     
    ui1=(1/c_j(1))*pij(:,1).*u_IMM(:,t); 
    ui2=(1/c_j(2))*pij(:,2).*u_IMM(:,t); 
    ui3=(1/c_j(3))*pij(:,3).*u_IMM(:,t);    %计算模型混合概率 
     
    % 计算各模型滤波初始化条件
    x01=x1_IMM*ui1(1)+x2_IMM*ui1(2)+x3_IMM*ui1(3); 
    x02=x1_IMM*ui2(1)+x2_IMM*ui2(2)+x3_IMM*ui2(3); 
    x03=x1_IMM*ui3(1)+x2_IMM*ui3(2)+x3_IMM*ui3(3);   %各模型滤波初始状态 
     
    P01=(P1_IMM+[x1_IMM-x01]*[x1_IMM-x01]')*ui1(1)+... 
        (P2_IMM+[x2_IMM-x01]*[x2_IMM-x01]')*ui1(2)+... 
        (P3_IMM+[x3_IMM-x01]*[x3_IMM-x01]')*ui1(3); 
    P02=(P1_IMM+[x1_IMM-x02]*[x1_IMM-x02]')*ui2(1)+... 
        (P2_IMM+[x2_IMM-x02]*[x2_IMM-x02]')*ui2(2)+... 
        (P3_IMM+[x3_IMM-x02]*[x3_IMM-x02]')*ui2(3); 
    P03=(P1_IMM+[x1_IMM-x03]*[x1_IMM-x03]')*ui3(1)+... 
        (P2_IMM+[x2_IMM-x03]*[x2_IMM-x03]')*ui3(2)+... 
        (P3_IMM+[x3_IMM-x03]*[x3_IMM-x03]')*ui3(3);  %各模型滤波初始状态协方差矩阵 
     
    %第二步--卡尔曼滤波 模型1卡尔曼滤波
    [x1_IMM,P1_IMM,r1_IMM,S1_IMM]=Kalman(x01,P01,z(:,t+1),F1,G,Q,H,R); 
    %模型2卡尔曼滤波
    [x2_IMM,P2_IMM,r2_IMM,S2_IMM]=Kalman(x02,P02,z(:,t+1),F2,G,Q,H,R); 
    %模型3卡尔曼滤波
    [x3_IMM,P3_IMM,r3_IMM,S3_IMM]=Kalman(x03,P03,z(:,t+1),F3,G,Q,H,R); 
     
    %第三步--模型概率更新
    [u_IMM(:,t+1)]=Model_P_up(r1_IMM,r2_IMM,r3_IMM,S1_IMM,S2_IMM,S3_IMM,c_j);      
    %第四步--模型综合
    
[x_pro_IMM(:,t+1),P_IMM(:,:,t+1)]=Model_mix(x1_IMM,x2_IMM,x3_IMM,P1_IMM,P2_IMM,P3_IMM,u_IMM(:,t)); 
end 
%=============================== 绘图 =============================== 目标轨迹
figure(1) 
plot(z_true(1,:),z_true(2,:)); grid on; hold on 
plot(x_pro_IMM(1,:),x_pro_IMM(3,:),'r'); 
plot(z(1,:),z(2,:),'*');  
hold off 
title('目标运动轨迹'); 
xlabel('x/m'); ylabel('y/m'); 
legend('真实值','滤波值','量测值'); 
text(z(1,1)+500,z(2,1),'t=1'); 
 
% 位置误差
figure(2) 
subplot(2,1,1); 
t=1:simTime; 
plot(t,abs(x_pro_IMM(1,t)-x(1,t)),'LineWidth',2);grid on 
title('x坐标位置跟踪误差'); 
xlabel('t/s'); ylabel('x-error/m'); 
 
subplot(2,1,2); 
t=1:simTime; 
plot(t,abs(x_pro_IMM(3,t)-x(3,t)),'LineWidth',2);grid on 
title('y坐标位置跟踪误差'); 
xlabel('t/s'); ylabel('y-error/m'); 
 
% 速度误差
figure(3) 
subplot(2,1,1); 
t=1:simTime; 
plot(t,abs(x_pro_IMM(2,t)-x(2,t)),'LineWidth',2);grid on 
title('x坐标速度跟踪误差'); 
xlabel('t/s'); ylabel('vx-error/m'); 
 
subplot(2,1,2); 
t=1:simTime; 
plot(t,abs(x_pro_IMM(4,t)-x(4,t)),'LineWidth',2);grid on 
title('y坐标速度跟踪误差'); 
xlabel('t/s'); ylabel('vy-error/m'); 
 
% 模型概率 figure(4)
plot(t,u_IMM(1,t),'k:',t,u_IMM(2,t),'r-.',t,u_IMM(3,t),'b--','LineWidth',2);grid on 
title('IMM算法模型概率曲线'); 
xlabel('t/s'); ylabel('模型概率'); 
legend('模型1','模型2','模型3'); 
end
 
function [u]=Model_P_up(r1,r2,r3,S1,S2,S3,c_j) 
%模型概率更新函数
 
%u  模型概率 r1 模型1预测误差 r2 模型2预测误差 r3 模型3预测误差 S1 模型1预测误差协方差矩阵 S2 模型2预测误差协方差矩阵
%S3 模型3预测误差协方差矩阵 c_j  模型混合概率
 
%计算似然函数
Lfun1=(1/sqrt(abs(2*pi*(det(S1)))))*exp((-1/2)*(r1'*inv(S1)*r1));   %Lfun1=1/(r1'*inv(S1)*r1); 
Lfun2=(1/sqrt(abs(2*pi*(det(S2)))))*exp((-1/2)*(r2'*inv(S2)*r2));   %Lfun2=1/(r2'*inv(S2)*r2); 
Lfun3=(1/sqrt(abs(2*pi*(det(S3)))))*exp((-1/2)*(r3'*inv(S3)*r3));    %Lfun3=1/(r3'*inv(S3)*r3); 
%归一化
Lfun11=Lfun1/(Lfun1+Lfun2+Lfun3); 
Lfun21=Lfun2/(Lfun1+Lfun2+Lfun3); 
Lfun31=Lfun3/(Lfun1+Lfun2+Lfun3); 
% 计算模型更新概率，即乘以上一时刻模型概率
c=[Lfun11,Lfun21,Lfun31]*c_j; 
% 再归一化
u=(1/c).*[Lfun11,Lfun21,Lfun31]'.*c_j; 
end

function [x_pro,P]=Model_mix(x1,x2,x3,P1,P2,P3,u) 
% 模型综合函数
%x_pro  状态综合值 P      综合协方差矩阵 x1     模型1状态估计值 x2     模型2状态估计值 x3
%模型3状态估计值 P1     模型1状态估计协方差矩阵 P2     模型2状态估计协方差矩阵 P3     模型3状态估计协方差矩阵 u
%模型转换概率
 % 按概率加权OK
x_pro=x1*u(1)+x2*u(2)+x3*u(3); 
P=(P1+[x1-x_pro]*[x1-x_pro]')*u(1)+... 
  (P2+[x2-x_pro]*[x2-x_pro]')*u(2)+... 
  (P3+[x3-x_pro]*[x3-x_pro]')*u(3); 
end
 
function [X,P,e,S]=Kalman(X_Forward,P_Forward,Z,A,G,Q,H,R) 
%卡尔曼滤波2012.2.27   IMM专用，参数略有不同 参数说明
%       Z--观测数据矢量 A--系统模型状态矩阵 G--系统模型噪声系数矩阵 Q--系统模型噪声方差 H--量测系数矩阵
%       R--量测模型噪声协方差 X_Forward--前次估计状态矢量 P_Forward--前次估计状态协方差矩阵
 
%       X--输出估计状态矢量 P--输出估计状态协方差矩阵 e--残差 S--残差协方差矩阵
 
% 预测
X_Pre=A*X_Forward; 
P_Pre=A*P_Forward*A'+G*Q*G'; 
 
% 增益矩阵
K=P_Pre*H'*inv(H*P_Pre*H'+R)'; 
 
% Pzz = H*P_Forward*H'+ R;                   %S(k+1/k+1) 新息协方差 Pxz =
% P_Forward*H' ;                       %状态与量测之间的协方差 K =
% P_Forward*H'*(inv(Pzz));               %K(k+1) 增益
     
e = Z-H*(A*X_Forward); 
S=H*P_Pre*H'+R;  %残差协方差矩阵 
 
% 修正滤波值和误差协方差阵
X=A*X_Forward+K*(Z-H*(A*X_Forward)); 
 
M=K*H; 
n=size(M); 
I=eye(n); 
P=(I-K*H)*P_Pre*(I-K*H)'+ K*R*K'; 
end
  