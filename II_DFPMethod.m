% author：zhangchunli
% University：DLUT
% Date：2022.04.22

x0 = [3,-1,0,1]';   % 初始点
maxError = 1e-4;    % 最大误差
f_test = fun([0,0,0,0]);
g_test = gfun([0,0,0,0]);
[x,val,k]=DFP(x0);
disp("最优解x的转置:x'")
disp(x')
disp("最优值：y")
disp(val)
disp("迭代次数：k")
disp(k)

function [x, val, k] = DFP(x0)
    k = 0;  % 迭代控制
    max_iteration = 20000;   % 最大迭代次数
    e = 1e-4;   % 允许的最大误差
    Hk=eye(4);  %初始Hk为对角阵
    while(k<max_iteration)
        gk=gfun(x0);%求梯度
        if(norm(gk)<=e)
            break;
        end  %判断是否满足停止迭代

        dk=-Hk*gk;%确定搜索方向，这个公式课本有DFP算法计算
        step = wolfepowell(x0,dk); %wolfepowell 线搜索
        x=x0+step*dk;
        g=gfun(x);%这个就是新的初始点的梯度
        sk=x-x0;%新的点减去原来的点
        yk=g-gk;%新的点的梯度减去旧梯度
        if(sk'*yk>0)
            Hk=Hk-(Hk*yk*yk'*Hk)/(yk'*Hk*yk)+(sk*sk')/(sk'*yk);%公式
        end
        x0=x;
        k=k+1;
    end
    x=x0;%break跳出后输出x
    val=fun(x);
end

function f=fun(x)   %x为n维向量
%需要求解最优值的目标函数
    f = (x(1)+10*x(2))^2+(sqrt(5)*(x(3)-x(4)))^2+((x(2)-2*x(3))^2)^2+(sqrt(10)*(x(1)-x(4))^2)^2;
end

function  g=gfun(x)
%目标函数的梯度
    g = [2*(x(1)+10*x(2))+40*(x(1)-x(4))^3, 20*(x(1)+10*x(2))+4*(x(2)-2*x(3))^3, 10*(x(3)-x(4))-8*(x(2)-2*x(3))^3, -10*(x(3)-x(4))-40*(x(1)-x(4))^3]';
end

function lambda = wolfepowell(Xk, sk)%定义lambda搜索方法
    c1=0.1;
    c2=0.5;
    lambda=1; 
    a=0; 
    b=Inf; 
    while (1)
        if ~(fun(Xk)-fun(Xk+lambda*sk)>=-c1*lambda*gfun(Xk)'*sk)
            b = lambda;
            lambda= (lambda+a)/2;
            continue
        end

        if ~(gfun(Xk+lambda*sk)'*sk >= c2*gfun(Xk)'*sk)
            a = lambda;
            lambda= min([2*lambda, (b+lambda)/2]);
            continue
        end
        break
    end
end



