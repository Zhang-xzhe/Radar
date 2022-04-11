%%%%      目标跟踪单元结课任务       %%%%



%  类似宏定义
rawx=200;                     %原始数据的x
rawy=200;                     %原始数据的y
rawz=40;                      %原始数据的帧数
basisnum=10;                  %提取每一帧信号点的编号
basist=40;                    %帧数
basisxy=2;                    %坐标
testgate=5;                   %是否为信号的门限
resultnum=3;                  %初定的轨迹编号
resultt=40;                   %帧数
resultxy=2;                   %坐标
mark=0;                       %检测标志位
tracknum=0;                   %记录总航迹数

%         第一，提取数据
basis_data=zeros(basisnum,basist,basisxy);            %存放信号数据
for c=1:rawz
    i=1;
    for a=1:rawx
        for b=1:rawy
            if  (raw_data(a,b,c)>testgate)   %信号的门限
            basis_data(i,c,1)=a;             %被判断为信号才可以存下来
            basis_data(i,c,2)=b;
            i=i+1;
            end
        end
    end
end



%          处理数据
dotnum=zeros(basisnum);               %确定每一帧有多少个信号点，把数目存起来
for i=1:basist
    for j=1:basisnum
        if  basis_data(j,i,1)~=0||basis_data(j,i,2)~=0
            dotnum(i)=dotnum(i)+1;
        end
    end
end

tracknum=dotnum(1);

%          确定航迹头
result_data=zeros(resultnum,resultt,resultxy);             %存储轨迹的数组
for j=1:dotnum(1)                                          %把第一帧的数据全部导入
     result_data(j,1,:)=basis_data(j,1,:);
end



%          航迹起始        取前n帧
%          1、规则法
startnum=5;                                   %确定取多少帧作为航迹起始
vmax=40;                                      %最大速度
t=1;                                          %相隔多久取一帧
next_aeraR=vmax*t;                            %运动半径
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


for c=1:3       %二维散点
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