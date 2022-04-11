% NNDA-FA ���ѡ��״����ݴ���Ӧ�á�P116
% ��ά�ռ�����ֱ���˶���״̬����ΪX=[x,vx,y,vy] x1=x0+vxT y1=y0+vyT
% ���棺 1���ı������������nc����ʽ��ȡ���ֶ����� 2���ı���������R=[r 0; 0 r]����r 3���ı��������λ��q��ƫ����ʵλ�õĳ̶�
% ���⣺ ÿ�β����Ӳ�ʱ��Sk���ȷ����
 
function nnda_test()
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
A_Model=[1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];     %����ģ�� 
H=[1 0 0 0;0 0 1 0];   %����ģ�� 
Q=1;              %ʵ�ʹ������� 
Q_Model=1;   %����ģ�͵Ĺ������� 
G = [T^2/2 0; T 0; 0 T^2/2; 0 T];   %������Ȩ���� 
r=20; 
R=[r 0; 0 r];   %�������� 
X0=[200;0;10000;-15];  %��ʼ״̬ 
X(:,1)=X0; 
Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
Zk(:,1)=H*X(:,1)+Vk; 
gama=16;            %�Ӳ����ò��� 
lamda=0.0004;     %��λ�������������� 
%************************************************
%          ��������
%************************************************
for i=2:1:simTime 
    X(:,i)=A*X(:,i-1);          % ��ʵ״̬ 
    Vk=[sqrt(r)*randn;sqrt(r)*randn]; 
    Zk(:,i)=H*X(:,i)+Vk;      %��������ֵ 
end 
%************************************************
%          NNSF��ʼ�� %************************************************
Xk_NNSF=[210;0;10100;-16];  %��ʼ״̬����ʵ��ֵ���в�� 
R11=r; R22=r; R12=0; R21=0; 
Pkk_NNSF=[R11 R11/T R12 R12/T; 
    R11/T 2*R11/T^2 R12/T 2*R12/T^2; 
    R21 R21/T R22 R22/T; 
    R21/T 2*R21/T^2 R22/T 2*R22/T^2];   %��ʼЭ���� 
Xkk =Xk_NNSF ; % X0 
Pkk = Pkk_NNSF ; 
X_Pre = Xk_NNSF ; 
P_Pre = Pkk_NNSF ; 
P=R; 
for i=1:1:simTime 
    %************************************************
    %          �����Ӳ�
    %************************************************
    Sk=H*P_Pre*H'+P;              %��ϢЭ���� 
    Av=pi*gama*sqrt(det(Sk));   % ����ȷ����������� ��*��*������SK 
    nc=floor(10*Av*lamda+1);    %�����Ӳ�������10Av*��+1 ȡ��   nc = 20 
    % disp(nc);       %��ʾ�Ӳ���Ŀ
    %�������
    q=sqrt(10*Av)/2;  %�м���� 
    q=q/10;  %��Ϊ��С��������ⲻ�ֲܷ�̫�㣬���������Ŀ�� 
    a=X(1,i)-q; 
    b=X(1,i)+q; 
    c=X(3,i)-q; 
    d=X(3,i)+q; 
    xi=a+(b-a)*rand(1,nc); 
    yi=c+(d-c)*rand(1,nc); 
    clear Z_Matrix;     %���ڴ������ 
    clear PZ_Matrix;   %���ڴ������ 
    for j=1:nc 
        Z_Matrix(:,j) = [xi(j);yi(j)];  %�Ӳ����⣺Z_Matrix�����ǰnc�� 
    end 
    Z_Matrix(:,nc+1)=Zk(:,i);  %��ʵ���⣺Z_Matrix����ĵ�nc+1�� 
     
    PZ_Matrix = cat(3);  %����������� 
    for j=1:1:nc 
        PZ_Matrix = cat(3,PZ_Matrix,[q,0; 0,q]);    %PZ_Matrixά����2*2*nc�� 
        %  �Ӳ������Ӧ�ķ��� PZ_Matrix(:,:,j)=[q,0; 0,q]  (j=1:nc)
    end 
    PZ_Matrix = cat(3,PZ_Matrix,R);  %�����ʹ�����Ӧ�ķ��PZ_Matrix(:,:,nc+1)=R 
    %************************************************
    %          NNDA����
    %************************************************
    Z_Predict = H*X_Pre;               %����Ԥ�� 
    PZ_Predict = H*P_Pre*H'+R ;   %��ϢЭ���� 
    [Z,P] = myNNDA(Z_Matrix, PZ_Matrix, Z_Predict, PZ_Predict);   % NNDA�����ع�������Ͷ�Ӧ����     
    Z_NNDA(:,i) = Z;  %����������洢 
    %************************************************
    %          �������˲�
    %************************************************
    [Xk,Pk,Kk]=Kalman(Xkk,Pkk,Z,A_Model,G,Q_Model,H,P); 
    Xkk=Xk; 
    Pkk=Pk; 
    % Ԥ��
    X_Pre=A_Model*Xkk; 
    P_Pre=A_Model*Pkk*A_Model'+G*Q_Model*G'; 
    %ȡ������״ֵ̬
    Ex_NNSF(i)=Xkk(1);     %x 
    Evx_NNSF(i)=Xkk(2);    %vx 
    Ey_NNSF(i)=Xkk(3);    %y 
    Evy_NNSF(i)=Xkk(4);    %vx 
    error1_NNSF(i)=Ex_NNSF(i)-X(1,i);%Pkk(1,1);    %error-x 
    error2_NNSF(i)=Ey_NNSF(i)-X(3,i);    %error-vx 
    error3_NNSF(i)=Evx_NNSF(i)-X(2,i);    %error-y 
    error4_NNSF(i)=Evy_NNSF(i)-X(4,i);    %error-vx 
end 
%************************************************
%          ��ͼ
%************************************************
i=1:simTime; 
%�켣
figure 
plot(X(1,i),X(3,i),'-','LineWidth',2); %��ʵֵ 
grid on; hold on 
plot(Ex_NNSF(1,i),Ey_NNSF(1,i),'r-','LineWidth',2); %�˲�ֵ 
plot(Zk(1,i),Zk(2,i),'*'); %ʵ�ʲ���ֵ 
plot(Z_NNDA(1,i),Z_NNDA(2,i),'o'); %�����ϵĲ���ֵ 
legend('��ʵֵ','�˲�ֵ','ʵ������','��������'); 
title('Ŀ���˶��켣'); xlabel('x/m'); ylabel('y/m'); 
text(X(1,1)+1,X(3,1)+5,'t=1'); 
 
%λ�����
figure 
subplot(211) 
plot(abs(error1_NNSF(i)),'LineWidth',2); grid on 
title('λ�����'); xlabel('t/s'); ylabel('error-x/m'); 
subplot(212) 
plot(abs(error3_NNSF(i)),'LineWidth',2); grid on 
xlabel('t/s'); ylabel('error-y/m'); 
 
%�ٶ����
figure 
subplot(211) 
plot(abs(error2_NNSF(i)),'LineWidth',2); grid on 
title('�ٶ����'); xlabel('t/s'); ylabel('error-vx/m/s'); subplot(212) 
plot(abs(error4_NNSF(i)),'LineWidth',2); grid on 
xlabel('t/s'); ylabel('error-vy/m/s'); 
end
 
function [Z,P] = myNNDA(Z_Matrix, PZ_Matrix, Z_Predict, PZ_Predict) 
% ���ڽ����ݹ������� ���룺 Z_Matrix�������ڵ���Ч����ֵ�������Ӳ�����ʵ���⣩ PZ_Matrix����Ч����ֵ��������
% Z_Predict������Ԥ��ֵ PZ_Predict������Ԥ��ֵ�������� ����� Z������ͳ�ƾ������ԭ������ϵ�����ֵ
% P�������ϵ�����ֵ��Ӧ��Э����
nm=size(Z_Matrix); 
n=nm(2);    % ��������Ч���������,������ 
for i=1:1:n 
    e(:,i)=Z_Matrix(:,i)-Z_Predict;  %ÿ��������Ԥ��ֵ�ľ��� 
    S(:,:,i)=PZ_Predict+PZ_Matrix(:,:,i);  %��ӦЭ���X��R��Q������������£� 
    D(:,i)=e(:,i)'*inv(S(:,:,i))*e(:,i); %ͳ�ƾ��� 
end 
Z=Z_Matrix(:,1) ; 
P=PZ_Matrix(:,:,1); 
d=D(:,1); 
index=1; 
for i=2:1:n 
    if D(:,i)<d 
        d=D(:,i); 
        Z=Z_Matrix(:,i) ; 
        P=PZ_Matrix(:,:,i); 
        index=i; 
    end 
end 
end

function [X,P,K]=Kalman(X_Forward,P_Forward,Z,A,G,Q,H,R) 
%�������˲� author:LW date:2011.06.08 ����˵��
%       Z--�۲�����ʸ��
 
%       A--ϵͳģ��״̬���� G--ϵͳģ������ϵ������ Q--ϵͳģ���������� H--����ϵ������ R--����ģ������Э����
%       X_Forward--ǰ�ι���״̬ʸ�� %       P_Forward--ǰ�ι���״̬Э�������
 
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
  