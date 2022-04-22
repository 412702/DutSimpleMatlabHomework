% author：zhangchunli
% University：DLUT
% Date：2022.04.22

%本人学号 22111091 ，按照作业要求设置参数a 和 G
%n = 291;
global n G b;
n = 291;
a = unidrnd(10, n, 1);
% G本身就是目标函数二阶梯度，即海森阵
G = a * a' + unidrnd(2) * eye(n);
b = 0.5 * G * ones(n, 1);

x0 = linspace(0, 0, n)';
[x,val,k]=dampNewton(x0);
disp("最优解x的转置:x'")
disp(x')
disp("最优值：y")
disp(val)
disp("迭代次数：k")
disp(k)

function [x, val, k] = dampNewton(x0)
    global G;
    k = 0;  % 迭代控制
    max_iteration = 1000;   % 最大迭代次数
    e = 1e-10;   % 允许的最大误差
    while(k<max_iteration)
        gk=gfun(x0);%求梯度
        if(norm(gk)<=e)
            break;
        end  %判断是否满足停止迭代

        dk=- G^-1 * gk;%确定搜索方向，阻尼牛顿法下降放下按公式来
        % 精确线搜索确定步长
        step = - (gk'*dk)/(dk'*G*dk);
        x=x0+step*dk;
        x0=x;
        k=k+1;
    end
    x=x0;%break跳出后输出x
    val=fun(x);
end

function f=fun(x)   %x为n维向量
%需要求解最优值的目标函数
    global G b;
    f = 0.5 * x' * G * x + b' * x;
end

function  g=gfun(x)
%目标函数的梯度
    global G b;
    g = G * x + b;
end

