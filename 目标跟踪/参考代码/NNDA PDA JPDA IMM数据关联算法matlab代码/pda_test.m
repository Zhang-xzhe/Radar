 % PDA-FA算法实现
% 何友《雷达数据处理及应用》P116
 
% 二维空间匀速直线运动，状态向量为X=[x,vx,y,vy] x1=x0+vxT y1=y0+vyT
 
% 仿真： 1、改变虚假量测数量nc：公式求取、手动设置 2、改变量测噪声R=[r 0; 0 r]，即r 3、改变虚假量测位置q，偏离真实位置的程度
% 4、关联概率计算
 
function pda_test()
clc; 
clear; 
close all; 
%************************************************
%          参数设置
%************************************************
I=eye(4); 
T = 1;                            %采样间隔 
simTime = 100 ;             %仿真步数 
A=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];   %实际模型：CV 
H=[1 0 0 0;0 0 1 0];   %测量模型 
Q=0;              %实际过程噪声 
G = [T^2/2 0; T 0; 0 T^2/2; 0 T];   %噪声加权矩阵 
r=200; 
R=[r 0; 0 r];   %量测噪声 
X0=[200;0;10000;-15];  %初始状态 
X(:,1)=X0; 
Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
Zk(:,1)=H*X(:,1)+Vk; 
gama=16; 
lamda=0.0004; 
%************************************************
%          量测生成
%************************************************
for i=2:1:simTime 
    X(:,i)=A*X(:,i-1);          % 真实状态 
    Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
    Zk(:,i)=H*X(:,i)+Vk;      %生成量测值 
end 
%************************************************
%          PDA初始化
%************************************************
Xk_PDA=[200;0;10100;-16];  %初始状态、与实际值略有差别 
R11=r; R22=r; R12=0; R21=0; 
Pkk_PDA=[R11 R11/T R12 R12/T;     R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
    R21 R21/T R22 R22/T; 
    R21/T 2*R21/T^2 R22/T 2*R22/T^2];   %初始协方差 
Xkk = Xk_PDA ; 
Pkk = Pkk_PDA; 
X_Pre = A*Xkk; 
P_Pre=A*Pkk*A'+G*Q*G'; 
P=R; 
for i=1:1:simTime 
    %************************************************
    %          产生杂波
    %************************************************
    % 量测确认区域面积
    Sk=H*P_Pre*H'+ P; 
    Av=pi*gama*sqrt(det(Sk)); 
    % 准备生成杂波数目
    nc=floor(10*Av*lamda+1);%设置杂波数量 
    q=sqrt(Av)/2;  %q=sqrt(10*Av)/2; 
    a=X(1,i)-q; 
    b=X(1,i)+q; 
    c=X(3,i)-q; 
    d=X(3,i)+q; 
    % 生成代表杂波的nc个虚假量测
    xi=a+(b-a)*rand(1,nc); 
    yi=c+(d-c)*rand(1,nc); 
    clear Z_Matrix; 
    clear PZ_Matrix; 
    for j=1:nc 
        Z_Matrix(:,j) = [xi(j);yi(j)]; 
    end 
    Z_Matrix(:,nc+1)=Zk(:,i); 
    PZ_Matrix = cat(3); 
    for j=1:1:nc 
        PZ_Matrix = cat(3,PZ_Matrix,[q,0;0,q]); 
    end 
    PZ_Matrix = cat(3,PZ_Matrix,R); 
    %************************************************
    %          PDA关联
    %************************************************
    Z_Predict = H*X_Pre; 
    PZ_Predict = H*P_Pre*H' ; 
    [Combine_Z,Combine_R]=PDA(Z_Matrix, PZ_Matrix, Z_Predict, PZ_Predict) ; % PDA 
    Z_PDA(:,i) = Combine_Z ; 
    %************************************************
    %          卡尔曼滤波
    %************************************************
    P=Combine_R; 
    [Xk_PDA,Pk_PDA,Kk_PDA]=Kalman(Xkk,Pkk,Combine_Z,A,G,Q,H,P); 
    Xkk=Xk_PDA;     Pkk=Pk_PDA; 
    % 预测
    X_Pre=A*Xkk; 
    P_Pre=A*Pkk*A'+G*Q*G'; 
    %出各个状态值
    Ex_PDA(i)=Xkk(1); 
    Evx_PDA(i)=Xkk(2); 
    Ey_PDA(i)=Xkk(3); 
    Evy_PDA(i)=Xkk(4); 
    error1_PDA(i)=Ex_PDA(i)-X(1,i);%Pkk(1,1); 
    error2_PDA(i)=Ey_PDA(i)-X(3,i);%Pkk(2,2); 
    error3_PDA(i)=Evx_PDA(i)-X(2,i);%Pkk(3,3); 
    error4_PDA(i)=Evy_PDA(i)-X(4,i);%Pkk(4,4); 
end 
%************************************************
%          绘图
%************************************************
i=1:simTime; 
figure 
plot(X(1,i),X(3,i),'-','LineWidth',2);    %真实值 
grid on; hold on 
plot(Ex_PDA(1,i),Ey_PDA(1,i),'r-','LineWidth',2);    %滤波值 
plot(Zk(1,i),Zk(2,i),'*');                    %实际测量值 
plot(Z_PDA(1,i),Z_PDA(2,i),'o');       %组合测量值 
legend('真实值','滤波值','实际量测','组合量测'); 
title('目标运动轨迹'); xlabel('x/m'); ylabel('y/m'); 
text(X(1,1)+1,X(3,1)+5,'t=1'); 
 
%位置误差
figure 
subplot(211) 
plot(abs(error1_PDA(i)),'LineWidth',2); grid on 
title('位置误差'); xlabel('t/s'); ylabel('error-x/m'); 
subplot(212) 
plot(abs(error3_PDA(i)),'LineWidth',2); grid on 
xlabel('t/s'); ylabel('error-y/m'); 
 
%速度误差
figure 
subplot(211) 
plot(abs(error2_PDA(i)),'LineWidth',2); grid on 
title('速度误差'); xlabel('t/s'); ylabel('error-vx/m/s'); 
subplot(212) 
plot(abs(error4_PDA(i)),'LineWidth',2); grid on 
xlabel('t/s'); ylabel('error-vy/m/s'); 
end
 
 
function [Combine_Z,Combine_R] = PDA(Z_Matrix, PZ_Matrix, Z_Predict, PZ_Predict) 
% 概率数据关联，杂波空间密度为泊松分布随机变量 输入： Z_Matrix：波门内的所有有效量测值 PZ_Matrix：有效量测值的误差方差阵
% Z_Predict：预测量测值 PZ_Predict：预测量测值的误差方差阵 输出： Combine_R为组合量测
% Combine_R：组合量测对应的协方差 中间变量： beta为正确关联概率
lamda=0.0004; 
Pd=1;               %检测概率，当不取1时，后面的a计算出来都是0 
Pg=0.9997;       %门限概率 
 
nm=size(Z_Matrix); 
n=nm(2);   % 量测数量 
m=nm(1);  % 测量维数 
 
for i=1:1:n 
    e(:,i)=Z_Matrix(:,i)-Z_Predict; 
    S(:,:,i)=PZ_Predict+PZ_Matrix(:,:,i);  %新息协方差 X、R、Q互不相关条件下 
    % 何友 计算方法P115 式（7.36） a(i)=exp((-1/2)*e(i)'*inv(S(i))*e(i));
    % bk(i)=lamda*sqrt((2*pi)*det(S(i)))*(1-Pd*Pg)/Pd; 杨万海P86式（3-5-7）
    a(i)=Pd*exp((-1/2)*(e(:,i)'*inv(S(:,:,i))*e(:,i)));          
    bk(i)=lamda*(sqrt(2*pi))^m*sqrt(det(S(:,:,i)))*(1-Pd);    
end 
for i=1:1:n 
    beta_i(i)=a(i)/(bk(i) + sum(a));     
end 
% 扩充正确关联概率，使得每一维量测都有对应的关联概率
beta = beta_i; 
for i=1:m-1 
    beta=[beta;beta_i]; 
end 
M = beta.*Z_Matrix; 
Combine_Z=sum(M',1); 
Combine_Z=Combine_Z'; 
Combine_R=0; 
for i=1:n 
    Combine_R = Combine_R + (beta(:,i)*beta(:,i)').*PZ_Matrix(:,:,i); 
end 
beta_i(n); 
end
  
function [X,P,K]=Kalman(X_Forward,P_Forward,Z,A,G,Q,H,R)
%卡尔曼滤波 2012.2.27 参数说明
%       Z--观测数据矢量
 
%       A--系统模型状态矩阵 G--系统模型噪声系数矩阵 Q--系统模型噪声方差 H--量测系数矩阵 R--量测模型噪声协方差
%       X_Forward--前次估计状态矢量 P_Forward--前次估计状态协方差矩阵
 
%       X--输出估计状态矢量 P--输出估计状态协方差矩阵
 
% 预测
X_Pre=A*X_Forward; 
P_Pre=A*P_Forward*A'+G*Q*G'; 
 
% 增益矩阵
K=P_Pre*H'*inv(H*P_Pre*H'+R)'; 
 
% Pzz = H*P_Forward*H'+ R;                   %S(k+1/k+1) 新息协方差 Pxz =
% P_Forward*H' ;                       %状态与量测之间的协方差 K =
% P_Forward*H'*(inv(Pzz));               %K(k+1) 增益
 
% 修正滤波值和误差协方差阵
X=A*X_Forward+K*(Z-H*(A*X_Forward)); 
 
M=K*H; 
n=size(M); 
I=eye(n); 
P=(I-K*H)*P_Pre*(I-K*H)'+ K*R*K'; 
end
 
 
   
 
 
 
 
  