%%%%      Ŀ����ٵ�Ԫ�������       %%%%



%  ���ƺ궨��
rawx=200;                     %ԭʼ���ݵ�x
rawy=200;                     %ԭʼ���ݵ�y
rawz=40;                      %ԭʼ���ݵ�֡��
basisnum=10;                  %��ȡÿһ֡�źŵ�ı��
basist=40;                    %֡��
basisxy=2;                    %����
testgate=5;                   %�Ƿ�Ϊ�źŵ�����
resultnum=3;                  %�����Ĺ켣���
resultt=40;                   %֡��
resultxy=2;                   %����
mark=0;                       %����־λ
tracknum=0;                   %��¼�ܺ�����

%         ��һ����ȡ����
basis_data=zeros(basisnum,basist,basisxy);            %����ź�����
for c=1:rawz
    i=1;
    for a=1:rawx
        for b=1:rawy
            if  (raw_data(a,b,c)>testgate)   %�źŵ�����
            basis_data(i,c,1)=a;             %���ж�Ϊ�źŲſ��Դ�����
            basis_data(i,c,2)=b;
            i=i+1;
            end
        end
    end
end



%          ��������
dotnum=zeros(basisnum);               %ȷ��ÿһ֡�ж��ٸ��źŵ㣬����Ŀ������
for i=1:basist
    for j=1:basisnum
        if  basis_data(j,i,1)~=0||basis_data(j,i,2)~=0
            dotnum(i)=dotnum(i)+1;
        end
    end
end

tracknum=dotnum(1);

%          ȷ������ͷ
result_data=zeros(resultnum,resultt,resultxy);             %�洢�켣������
for j=1:dotnum(1)                                          %�ѵ�һ֡������ȫ������
     result_data(j,1,:)=basis_data(j,1,:);
end



%          ������ʼ        ȡǰn֡
%          1������
startnum=5;                                   %ȷ��ȡ����֡��Ϊ������ʼ
vmax=40;                                      %����ٶ�
t=1;                                          %������ȡһ֡
next_aeraR=vmax*t;                            %�˶��뾶
for a=1:startnum
        for c=1:dotnum(a+1)
            x1=result_data(c,a+1,1);
            y1=result_data(c,a+1,2);
           for b=1:tracknum
                x=result_data(b,a,1);
                y=result_data(b,a,2);
                if  (x1-x)^2+(y1-y)^2<vmax^2
                    result_data(b,a+1,1)=x1;
                    result_data(b,a+1,2)=y1;
                    mark=1;
                end
           end
            if mark==0
                tracknum=tracknum+1;
                result_data(tracknum,a+1,1)=x1;
                result_data(tarcknum,a+1,2)=y1;
            end
            mark=0;
        end
end


for c=1:3       %��άɢ��
   % figure;
    for a=1:40
        
            
           
            plot(result_data(c,a,1),result_data(c,a,2),'o');
            xlim([0,200]);
            ylim([0,200]);
            pause(0.5)
            hold on
            
        
    end
    % saveas(gcf,picname);
end 