%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Probabilistic Data Association Filter    http://wenku.baidu.com/link?url=Q7B63MhiQWiUfRUMZLTHSUPjqp6DBpA1msEbxd9HuJtR5njHEfFY7q5i2n1rn0e-ofvbo9Bt5f-_fjTQH6OJgPcU33stOE-ex-BtAxwJcyO
%Writed by Liangqun Li  
%Date:2006.4.21 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
tic 
n=50;
Pd=1;                                                    %检测概率 
Pg=0.99;                                                 %正确量测落入跟踪门内得概率 
g_sigma=9.21;                                            %门限 
gamma=1;                                                 %每一个单位面积内产生一个杂波 
Target_measurement=zeros(2,n);                           %目标观测互联存储矩阵 
target_delta=0.15;                                       %目标对应的观测标准差    
data_measurement=zeros(2,n);                             %观测存储矩阵,n采样次数    行---代表目标 列---代表传感器                       
P=zeros(4,4);                                            %滤波时协方差更新 
x_filter=zeros(4,n);                                     %存储所有六个目标的各时刻的滤波值 
A = [1 T 0 0;0 1 0 0;0 0 1 T;0 0 0 1];                   %状态转移矩阵 
C = [1 0 0 0;0 0 1 0];  
R=[target_delta^2  0;0  target_delta^2];                 % 
G=[T^2/2 0;0 T^2/2;T 0;0 T]; 
P=[target_delta^2 0 0 0;0 0.01 0 0;0 0 target_delta^2 0;0 0 0 0.01];%初始化协方差 
x_filter1=zeros(4,n,MC_number);                          %MC_number次montle carlo仿真所得全部结果存储矩阵 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%主程序 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
for M=1:MC_number 
%产生路径 
Noise=[]; 
for i=1:n 
            data_measurement(1,i)=data_measurement1(1,i)+randn(1)*target_delta; 
            data_measurement(2,i)=data_measurement1(2,i)+randn(1)*target_delta;                 %传感器观测的位置  
end 
NOISE=[]; 
%滤波开始 
for t=1:n 
%     Noise=[]; 
%     y=[]; 
%     for i=1:c 
%         Noise_x=V(i,1)+target_delta(i)*3.5-2*rand(1,mm)*target_delta(i)*3.5; 
%         Noise_y=V(i,2)+target_delta(i)*3.5-2*rand(1,mm)*target_delta(i)*3.5;      %产生杂波 
%         Noise1=[Noise_x ;Noise_y]; 
%         Noise=[Noise Noise1]; 
%     end 
%     b=zeros(1,2); 
%     b(1)=data_measurement(1,1,t); 
%     b(2)=data_measurement(1,2,t); 
%     y=[Noise b']; 
%     b(1)=data_measurement(2,1,t); 
%     b(2)=data_measurement(2,2,t); 
%     y=[y b'];                                                                     %杂波、观测都放在y中 
 
    if t~=1 
        x_predic = A*x_filter(:,t-1);                                  %用前一时刻的滤波值来预测当前的值  
    else 
        x_predic = target_position';                                   %第一次采样我们用真实位置当预测值  
    end 
    P_predic = A*P*A'+G*Q*G'; 
    Z_predic = C*x_predic; 
    S = C*P_predic*C'+ R; 
    K = P_predic*C'*inv(S);                                            %增益 
    ellipse_Volume=pi*g_sigma*sqrt(det(S));                            %计算椭球体积，这里算的是面积    
    number_returns=floor(10*ellipse_Volume*gamma+1);                   %错误回波数 
    side=sqrt((10*ellipse_Volume*gamma+1)/gamma)/2;                    %求出正方行边长的二分之一 
    Noise_x=x_predic(1)+side-2*rand(1,number_returns)*side;            %在预测值周围产生多余回波 
    Noise_y=x_predic(3)+side-2*rand(1,number_returns)*side; 
    Noise=[Noise_x ;Noise_y]; 
    NOISE=[NOISE Noise]; 
    % 
    b=zeros(1,2); 
    b(1)=data_measurement(1,t); 
    b(2)=data_measurement(2,t); 
    y1=[Noise b'];                                                      %将接收到的所有的回波存在y1中 
    y=[]; 
    d=[]; 
    m=0; 
    for j=1:(number_returns+1) 
        d=y1(:,j)-Z_predic; 
        D=d'*inv(S)*d; 
        if D<=g_sigma 
            y=[y y1(:,j)];                                              %把落入跟踪门中的所有回波放入y中 
            m=m+1;                                                      %记录观测个数 
        end 
    end 
     
    %m=0表示无有效回波 
    Bk=gamma*2*pi*sqrt(det(S))*(1-Pd*Pg)/Pd;                            %算b0       
    if m==0 
       x_filter(:,t)= x_predic; 
       P=P_predic;                                                      %无回波的情况 
   else         
    E=zeros(1,m); 
    belta=zeros(1,m); 
    for i=1:m 
        a=(y(:,i)-Z_predic)'*inv(S)*(y(:,i)-Z_predic); 
        E(i)=E(i)+exp(-a/2); 
    end 
    belta0=Bk/(Bk+sum(E));                                                                                                         %无回波时的关联概率 
    v=zeros(2,1); 
    v1=zeros(2,2); 
    for i=1:m 
        belta(i)=E(i)/(Bk+sum(E));                                          %算关联概率 
        v=v+belta(i)*(y(:,i)-Z_predic); 
        v1=v1+belta(i)*(y(:,i)-Z_predic)*(y(:,i)-Z_predic)'; 
    end 
    x_filter(:,t)= x_predic + K*v; 
%算协方差 
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
% %画图 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%R=sum(U1,4)/MC_number;   
x_filter=sum(x_filter1,3)/MC_number;                                %滤波值作平均 
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
        a(1,j)=sqrt((x_filter(1,j)-data_measurement1(1,j))^2+(x_filter(3,j)-data_measurement1(2,j))^2);%最小均方误差 
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
%     a(1,j)=sqrt((x_filter(1,1,j)-data_measurement1(1,1,j))^2+(x_filter(3,1,j)-data_measurement1(1,2,j))^2);%最小均方误差 
%     c1(1,j)=c1(1,j)+a(1,j); 
% end 
% % axis([0 20 0 11 ]) 
% xlabel('times'),ylabel('Mean Estimated Errors/m'); 
% legend('target1','target2') 
% % 
% c1=zeros(1,n); 
% for j=1:n 
%     a(1,j)=sqrt((x_filter(2,1,j)-target_position(2,1))^2+(x_filter(4,1,j)-target_position(4,1))^2);%最小均方误差 
%     c1(1,j)=c1(1,j)+a(1,j); 
% end 
% figure(4) 
% plot(1:n,c1(1,:),'r:')  
% hold on 
% c1=zeros(1,n); 
% for j=1:n 
%     a(1,j)=sqrt((x_filter(2,2,j)-target_position(2,2))^2+(x_filter(4,2,j)-target_position(4,2))^2);%最小均方误差 
%     c1(1,j)=c1(1,j)+a(1,j); 
% plot(1:n,c1(1,:),'r-')  
% xlabel('times'),ylabel('Mean Estimated Errors m/s'); 
% legend('target1','target2') 