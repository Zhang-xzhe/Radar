%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Probabilistic Data Association Filter    http://wenku.baidu.com/link?url=Q7B63MhiQWiUfRUMZLTHSUPjqp6DBpA1msEbxd9HuJtR5njHEfFY7q5i2n1rn0e-ofvbo9Bt5f-_fjTQH6OJgPcU33stOE-ex-BtAxwJcyO
%Writed by Liangqun Li  
%Date:2006.4.21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
tic 
n=50;
Pd=1;                                                    %������ 
Pg=0.99;                                                 %��ȷ��������������ڵø��� 
g_sigma=9.21;                                            %���� 
gamma=1;                                                 %ÿһ����λ����ڲ���һ���Ӳ� 
Target_measurement=zeros(2,n);                           %Ŀ��۲⻥���洢���� 
target_delta=0.15;                                       %Ŀ���Ӧ�Ĺ۲��׼��    
data_measurement=zeros(2,n);                             %�۲�洢����,n��������    ��---����Ŀ�� ��---��������                       
P=zeros(4,4);                                            %�˲�ʱЭ������� 
x_filter=zeros(4,n);                                     %�洢��������Ŀ��ĸ�ʱ�̵��˲�ֵ 
A = [1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];                   %״̬ת�ƾ��� 
C = [1 0 0 0;0 0 1 0];  
R=[target_delta^2  0;0  target_delta^2];                 % 
G=[T^2/2 0;0 T^2/2;T 0;0 T]; 
P=[target_delta^2 0 0 0;0 0.01 0 0;0 0 target_delta^2 0;0 0 0 0.01];%��ʼ��Э���� 
x_filter1=zeros(4,n,MC_number);                          %MC_number��montle carlo��������ȫ������洢���� 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%������ 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
for M=1:MC_number 
%����·�� 
Noise=[]; 
for i=1:n 
            data_measurement(1,i)=data_measurement1(1,i)+randn(1)*target_delta; 
            data_measurement(2,i)=data_measurement1(2,i)+randn(1)*target_delta;                 %�������۲��λ��  
end 
NOISE=[]; 
%�˲���ʼ 
for t=1:n 
%     Noise=[]; 
%     y=[]; 
%     for i=1:c 
%         Noise_x=V(i,1)+target_delta(i)*3.5-2*rand(1,mm)*target_delta(i)*3.5; 
%         Noise_y=V(i,2)+target_delta(i)*3.5-2*rand(1,mm)*target_delta(i)*3.5;      %�����Ӳ� 
%         Noise1=[Noise_x ;Noise_y]; 
%         Noise=[Noise Noise1]; 
%     end 
%     b=zeros(1,2); 
%     b(1)=data_measurement(1,1,t); 
%     b(2)=data_measurement(1,2,t); 
%     y=[Noise b']; 
%     b(1)=data_measurement(2,1,t); 
%     b(2)=data_measurement(2,2,t); 
%     y=[y b'];                                                                     %�Ӳ����۲ⶼ����y�� 
 
    if t~=1 
        x_predic = A*x_filter(:,t-1);                                  %��ǰһʱ�̵��˲�ֵ��Ԥ�⵱ǰ��ֵ  
    else 
        x_predic = target_position';                                   %��һ�β�����������ʵλ�õ�Ԥ��ֵ  
    end 
    P_predic = A*P*A'+G*Q*G'; 
    Z_predic = C*x_predic; 
    S = C*P_predic*C'+ R; 
    K = P_predic*C'*inv(S);                                            %���� 
    ellipse_Volume=pi*g_sigma*sqrt(det(S));                            %�������������������������    
    number_returns=floor(10*ellipse_Volume*gamma+1);                   %����ز��� 
    side=sqrt((10*ellipse_Volume*gamma+1)/gamma)/2;                    %��������б߳��Ķ���֮һ 
    Noise_x=x_predic(1)+side-2*rand(1,number_returns)*side;            %��Ԥ��ֵ��Χ��������ز� 
    Noise_y=x_predic(3)+side-2*rand(1,number_returns)*side; 
    Noise=[Noise_x ;Noise_y]; 
    NOISE=[NOISE Noise]; 
    % 
    b=zeros(1,2); 
    b(1)=data_measurement(1,t); 
    b(2)=data_measurement(2,t); 
    y1=[Noise b'];                                                      %�����յ������еĻز�����y1�� 
    y=[]; 
    d=[]; 
    m=0; 
    for j=1:(number_returns+1) 
        d=y1(:,j)-Z_predic; 
        D=d'*inv(S)*d; 
        if D<=g_sigma 
            y=[y y1(:,j)];                                              %������������е����лز�����y�� 
            m=m+1;                                                      %��¼�۲���� 
        end 
    end 
     
    %m=0��ʾ����Ч�ز� 
    Bk=gamma*2*pi*sqrt(det(S))*(1-Pd*Pg)/Pd;                            %��b0       
    if m==0 
       x_filter(:,t)= x_predic; 
       P=P_predic;                                                      %�޻ز������ 
   else         
    E=zeros(1,m); 
    belta=zeros(1,m); 
    for i=1:m 
        a=(y(:,i)-Z_predic)'*inv(S)*(y(:,i)-Z_predic); 
        E(i)=E(i)+exp(-a/2); 
    end 
    belta0=Bk/(Bk+sum(E));                                                                                                         %�޻ز�ʱ�Ĺ������� 
    v=zeros(2,1); 
    v1=zeros(2,2); 
    for i=1:m 
        belta(i)=E(i)/(Bk+sum(E));                                          %��������� 
        v=v+belta(i)*(y(:,i)-Z_predic); 
        v1=v1+belta(i)*(y(:,i)-Z_predic)*(y(:,i)-Z_predic)'; 
    end 
    x_filter(:,t)= x_predic + K*v; 
%��Э���� 
    Pc=(eye(4)-K*C)*P_predic; 
    PP=K*(v1-v*v')*K'; 
    P=belta0*P_predic+(1-belta0)*Pc+PP; 
    end 
%     Z=y(:,num); 
%     K = P_predic*C'*inv(S); 
%     x_filter(:,i,t)= x_predic + K*(Z - Z_predic); 
%     P(:,:,i)= P_predic - K*S*K'; 
%     x_predic = A*x_filter(:,i,t); 
%     V(i,1)=x_predic(1);V(i,2)=x_predic(3); 
    x_filter1(:,t,M)=x_filter(:,t); 
% end% end 
end 
end 
toc 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% %��ͼ 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%R=sum(U1,4)/MC_number;   
x_filter=sum(x_filter1,3)/MC_number;                                %�˲�ֵ��ƽ�� 
plot(x_filter(1,:),x_filter(3,:),'r+'),hold on 
plot(data_measurement(1,:),data_measurement(2,:),'*') 
plot(data_measurement1(1,:),data_measurement1(2,:),'-') 
% plot(NOISE(1,:),NOISE(2,:),'.') 
axis([0 30 1 7]) 
a=zeros(1,n); 
xlabel('x(km)'),ylabel('y(km)'); 
legend('estimated position','measurements','target position',4) 
figure 
for j=1:n 
        a(1,j)=sqrt((x_filter(1,j)-data_measurement1(1,j))^2+(x_filter(3,j)-data_measurement1(2,j))^2);%��С������� 
end 
plot(1:n,a(1,:),'r:')  
xlabel('t(s)'),ylabel('estimated errors(km)'); 
 
% a=zeros(1,n); 
% b=zeros(1,n); 
% figure(1) 
% for i=1:c 
%     a=zeros(1,n); 
%     b=zeros(1,n); 
%     for j=1:n 
%         a(j)=data_measurement1(i,1,j);b(j)=data_measurement1(i,2,j); 
%     end 
%     plot(a(:),b(:),'-'),hold on 
% end 
% % axis([5900 8500 5900 6200 ]) 
% for i=1:c 
 
 
%     for j=1:n 
%         a(j)=x_filter(1,i,j);b(j)=x_filter(3,i,j); 
%     end 
% if i==1 
%     plot(a(:),b(:),'m:') 
% else  
%     plot(a(:),b(:),':') 
% end 
% end 
% xlabel('x/m'),ylabel('y/m'); 
% axis([5900 9100 6150 6450 ]) 
% a=zeros(c,n);b=zeros(c,n);c1=zeros(1,n);MEASURE=zeros(2,c,n); 
% % for j=1:n 
% %     for i=1:c 
% %         hh=0;hh1=0; 
% %         a(i,j)=sqrt((x_filter(1,i,j)-data_measurement1(i,1,j))^2+(x_filter(3,i,j)-data_measurement1(i,2,j))^2); 
% %         for k=1:5 
% %             hh=data_measurement(i,1,k,j)+hh;hh1=data_measurement(i,2,k,j)+hh1; 
% %             MEASURE(:,i,j)=[hh/5 hh1/5]'; 
% %         end 
% %      b(i,j)=sqrt((MEASURE(1,i,j)-data_measurement1(i,1,j))^2+(MEASURE(2,i,j)-data_measurement1(i,2,j))^2); 
% %         end 
% %     c1(1,j)=sum(a(i,j))/sum(b(i,j)); 
% % end 
% for j=1:n 
%     a(1,j)=sqrt((x_filter(1,1,j)-data_measurement1(1,1,j))^2+(x_filter(3,1,j)-data_measurement1(1,2,j))^2);%��С������� 
%     c1(1,j)=c1(1,j)+a(1,j); 
% end 
% % axis([0 20 0 11 ]) 
% xlabel('times'),ylabel('Mean Estimated Errors/m'); 
% legend('target1','target2') 
% % 
% c1=zeros(1,n); 
% for j=1:n 
%     a(1,j)=sqrt((x_filter(2,1,j)-target_position(2,1))^2+(x_filter(4,1,j)-target_position(4,1))^2);%��С������� 
%     c1(1,j)=c1(1,j)+a(1,j); 
% end 
% figure(4) 
% plot(1:n,c1(1,:),'r:')  
% hold on 
% c1=zeros(1,n); 
% for j=1:n 
%     a(1,j)=sqrt((x_filter(2,2,j)-target_position(2,2))^2+(x_filter(4,2,j)-target_position(4,2))^2);%��С������� 
%     c1(1,j)=c1(1,j)+a(1,j); 
% plot(1:n,c1(1,:),'r-')  
% xlabel('times'),ylabel('Mean Estimated Errors m/s'); 
% legend('target1','target2') 