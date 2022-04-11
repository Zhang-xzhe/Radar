a=0;
sigma=1;
sample=100000;
r = mvnrnd(a,sigma,sample);
subplot(2,1,1);
plot(r);
x = linspace(min(r), max(r), 1000); 
y = hist(r,x); 
y = y/sample; 
subplot(2,1,2);
bar(x,y);
hold on;
s = 0;
for i=2:length(x) 
    s=[s,trapz(x(1:i),y(1:i))];
end
%     figure;
%     plot(x,s,x,s,'*')
%     a= polyfit(x,y,2);
%     yy= polyval(a, x);
% %     plot(x, yy); %拟合后的曲线


 
