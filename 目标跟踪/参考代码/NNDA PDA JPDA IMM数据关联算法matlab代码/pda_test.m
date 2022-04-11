 % PDA-FA�㷨ʵ��
% ���ѡ��״����ݴ���Ӧ�á�P116
 
% ��ά�ռ�����ֱ���˶���״̬����ΪX=[x,vx,y,vy] x1=x0+vxT y1=y0+vyT
 
% ���棺 1���ı������������nc����ʽ��ȡ���ֶ����� 2���ı���������R=[r 0; 0 r]����r 3���ı��������λ��q��ƫ����ʵλ�õĳ̶�
% 4���������ʼ���
 
function pda_test()
clc; 
clear; 
close all; 
%************************************************
%          ��������
%************************************************
I=eye(4); 
T = 1;                            %������� 
simTime = 100 ;             %���沽�� 
A=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];   %ʵ��ģ�ͣ�CV 
H=[1 0 0 0;0 0 1 0];   %����ģ�� 
Q=0;              %ʵ�ʹ������� 
G = [T^2/2 0; T 0; 0 T^2/2; 0 T];   %������Ȩ���� 
r=200; 
R=[r 0; 0 r];   %�������� 
X0=[200;0;10000;-15];  %��ʼ״̬ 
X(:,1)=X0; 
Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
Zk(:,1)=H*X(:,1)+Vk; 
gama=16; 
lamda=0.0004; 
%************************************************
%          ��������
%************************************************
for i=2:1:simTime 
    X(:,i)=A*X(:,i-1);          % ��ʵ״̬ 
    Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
    Zk(:,i)=H*X(:,i)+Vk;      %��������ֵ 
end 
%************************************************
%          PDA��ʼ��
%************************************************
Xk_PDA=[200;0;10100;-16];  %��ʼ״̬����ʵ��ֵ���в�� 
R11=r; R22=r; R12=0; R21=0; 
Pkk_PDA=[R11 R11/T R12 R12/T;     R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
    R21 R21/T R22 R22/T; 
    R21/T 2*R21/T^2 R22/T 2*R22/T^2];   %��ʼЭ���� 
Xkk = Xk_PDA ; 
Pkk = Pkk_PDA; 
X_Pre = A*Xkk; 
P_Pre=A*Pkk*A'+G*Q*G'; 
P=R; 
for i=1:1:simTime 
    %************************************************
    %          �����Ӳ�
    %************************************************
    % ����ȷ���������
    Sk=H*P_Pre*H'+ P; 
    Av=pi*gama*sqrt(det(Sk)); 
    % ׼�������Ӳ���Ŀ
    nc=floor(10*Av*lamda+1);%�����Ӳ����� 
    q=sqrt(Av)/2;  %q=sqrt(10*Av)/2; 
    a=X(1,i)-q; 
    b=X(1,i)+q; 
    c=X(3,i)-q; 
    d=X(3,i)+q; 
    % ���ɴ����Ӳ���nc���������
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
    %          PDA����
    %************************************************
    Z_Predict = H*X_Pre; 
    PZ_Predict = H*P_Pre*H' ; 
    [Combine_Z,Combine_R]=PDA(Z_Matrix, PZ_Matrix, Z_Predict, PZ_Predict) ; % PDA 
    Z_PDA(:,i) = Combine_Z ; 
    %************************************************
    %          �������˲�
    %************************************************
    P=Combine_R; 
    [Xk_PDA,Pk_PDA,Kk_PDA]=Kalman(Xkk,Pkk,Combine_Z,A,G,Q,H,P); 
    Xkk=Xk_PDA;     Pkk=Pk_PDA; 
    % Ԥ��
    X_Pre=A*Xkk; 
    P_Pre=A*Pkk*A'+G*Q*G'; 
    %������״ֵ̬
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
%          ��ͼ
%************************************************
i=1:simTime; 
figure 
plot(X(1,i),X(3,i),'-','LineWidth',2);    %��ʵֵ 
grid on; hold on 
plot(Ex_PDA(1,i),Ey_PDA(1,i),'r-','LineWidth',2);    %�˲�ֵ 
plot(Zk(1,i),Zk(2,i),'*');                    %ʵ�ʲ���ֵ 
plot(Z_PDA(1,i),Z_PDA(2,i),'o');       %��ϲ���ֵ 
legend('��ʵֵ','�˲�ֵ','ʵ������','�������'); 
title('Ŀ���˶��켣'); xlabel('x/m'); ylabel('y/m'); 
text(X(1,1)+1,X(3,1)+5,'t=1'); 
 
%λ�����
figure 
subplot(211) 
plot(abs(error1_PDA(i)),'LineWidth',2); grid on 
title('λ�����'); xlabel('t/s'); ylabel('error-x/m'); 
subplot(212) 
plot(abs(error3_PDA(i)),'LineWidth',2); grid on 
xlabel('t/s'); ylabel('error-y/m'); 
 
%�ٶ����
figure 
subplot(211) 
plot(abs(error2_PDA(i)),'LineWidth',2); grid on 
title('�ٶ����'); xlabel('t/s'); ylabel('error-vx/m/s'); 
subplot(212) 
plot(abs(error4_PDA(i)),'LineWidth',2); grid on 
xlabel('t/s'); ylabel('error-vy/m/s'); 
end
 
 
function [Combine_Z,Combine_R] = PDA(Z_Matrix, PZ_Matrix, Z_Predict, PZ_Predict) 
% �������ݹ������Ӳ��ռ��ܶ�Ϊ���ɷֲ�������� ���룺 Z_Matrix�������ڵ�������Ч����ֵ PZ_Matrix����Ч����ֵ��������
% Z_Predict��Ԥ������ֵ PZ_Predict��Ԥ������ֵ�������� ����� Combine_RΪ�������
% Combine_R����������Ӧ��Э���� �м������ betaΪ��ȷ��������
lamda=0.0004; 
Pd=1;               %�����ʣ�����ȡ1ʱ�������a�����������0 
Pg=0.9997;       %���޸��� 
 
nm=size(Z_Matrix); 
n=nm(2);   % �������� 
m=nm(1);  % ����ά�� 
 
for i=1:1:n 
    e(:,i)=Z_Matrix(:,i)-Z_Predict; 
    S(:,:,i)=PZ_Predict+PZ_Matrix(:,:,i);  %��ϢЭ���� X��R��Q������������� 
    % ���� ���㷽��P115 ʽ��7.36�� a(i)=exp((-1/2)*e(i)'*inv(S(i))*e(i));
    % bk(i)=lamda*sqrt((2*pi)*det(S(i)))*(1-Pd*Pg)/Pd; ����P86ʽ��3-5-7��
    a(i)=Pd*exp((-1/2)*(e(:,i)'*inv(S(:,:,i))*e(:,i)));          
    bk(i)=lamda*(sqrt(2*pi))^m*sqrt(det(S(:,:,i)))*(1-Pd);    
end 
for i=1:1:n 
    beta_i(i)=a(i)/(bk(i) + sum(a));     
end 
% ������ȷ�������ʣ�ʹ��ÿһά���ⶼ�ж�Ӧ�Ĺ�������
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
%�������˲� 2012.2.27 ����˵��
%       Z--�۲�����ʸ��
 
%       A--ϵͳģ��״̬���� G--ϵͳģ������ϵ������ Q--ϵͳģ���������� H--����ϵ������ R--����ģ������Э����
%       X_Forward--ǰ�ι���״̬ʸ�� P_Forward--ǰ�ι���״̬Э�������
 
%       X--�������״̬ʸ�� P--�������״̬Э�������
 
% Ԥ��
X_Pre=A*X_Forward; 
P_Pre=A*P_Forward*A'+G*Q*G'; 
 
% �������
K=P_Pre*H'*inv(H*P_Pre*H'+R)'; 
 
% Pzz = H*P_Forward*H'+ R;                   %S(k+1/k+1) ��ϢЭ���� Pxz =
% P_Forward*H' ;                       %״̬������֮���Э���� K =
% P_Forward*H'*(inv(Pzz));               %K(k+1) ����
 
% �����˲�ֵ�����Э������
X=A*X_Forward+K*(Z-H*(A*X_Forward)); 
 
M=K*H; 
n=size(M); 
I=eye(n); 
P=(I-K*H)*P_Pre*(I-K*H)'+ K*R*K'; 
end
 
 
   
 
 
 
 
  