%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%generate data of target trajectory 
%Writed by Liangqun Li  
%Date:2006.4.21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear all 
n=50;                                                    %�������� 
T=1;                                                     %TΪ������� 
MC_number=10;                                            %monte carlos run times 
target_position=[1.5 0.5 1.5  0.1];                      %Ŀ�����ʼλ�ú��ٶ�                    
data_measurement1=ze
ros(2,n);                            %data_measurement�۲�ֵ����,data_measurement1ʵ��λ�þ���     
Q=[0.0004 0;0 0.0004]; 
Qdelta=sqrt(Q(1,1)); 
data_measurement1(:,1,1)=target_position(1); 
data_measurement1(:,2,1)=target_position(3); 
for i=2:n 
        if i~=1 
            data_measurement1(1,i)=data_measurement1(1,1)+T*(i-1)*target_position(2)+rand(1)*Qdelta;            
            data_measurement1(2,i)=data_measurement1(2,1)+T*(i-1)*target_position(4)+rand(1)*Qdelta;   %ʵ��λ�� �������ٶ� 
        end 
end 
plot(data_measurement1(1,:),data_measurement1(2,:),'-'); 
axis([0 30 1 7])