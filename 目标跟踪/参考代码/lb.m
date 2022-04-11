function [k,b,vx,vy] = lb(X, Y, num, T)
%   y=kx+b
    x_ave=0;
    y_ave=0;
    fenzi=0;
    fenmu=0;
    if(num>0)
        for i=1:num
            x_ave=x_ave+X(i);
            y_ave=y_ave+Y(i);
        end
        x_ave=x_ave/num;
        y_ave=y_ave/num;
        for i=1:num
            fenzi=fenzi+(X(i)-x_ave)*(Y(i)-y_ave);
            fenmu=fenmu+(X(i)-x_ave)*(X(i)-x_ave);
        end
        k=fenzi/fenmu;
        b=y_ave-x_ave*k;
        
%         dx=X(num)-X(1);
%         vx=dx/(num*T);
%         vy=k*dx/(num*T);

        if(mod(num,2)==0)
            n=num/2;
        else
            n=(num-1)/2;
        end
        addx=0;
        addy=0;
        for i=1:(num-n)
            addx=addx+X(i+n)-X(i);
            addy=addy+Y(i+n)-Y(i);
        end
        vx=addx/(n*(num-n));
        vy=addy/(n*(num-n));
    end
end