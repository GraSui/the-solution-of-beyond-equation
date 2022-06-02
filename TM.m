%Set x=w_p/(c*k_gx)
syms x y yy
P=[]
options=optimoptions("fsolve","Algorithm","trust-region")
f=str2func(['@(x,y)',vectorize(F(x,y))])
y=0.0005:0.002:1.6005;Y={};
for x=0:0.02:1.2
    %index_y=y>x   %select the conducting modes
    %y=y(index_y)
    g=str2func(['@(yy)',vectorize(f(x,yy))])
    f_scatter=g(y)
    [g_scatter,index]=sort(abs(f_scatter))
    y_index=y(index(1:55))
    y_index=sort(y_index)
    c=1;C={};count=1
    for i=2:length(y_index)
        if y_index(i)-y_index(i-1)>=0.0021
            C=[C,y_index(c:i-1)]
            c=i
            count=count+1
        end
    end
    C=[C,y_index(c:length(y_index))]  
    % selection for the second time
    Y_index2=[]
    for i=1:count
        D=C{i}
        f_scatter2=g(D)
        [g_scatter2,index2]=sort(abs(f_scatter2))
        y_index2=D(index2(1))
        Y_index2=[Y_index2,y_index2]
    end

    y_x=fsolve(g,Y_index2,options)
    index_y=abs(y_x-x)>0.0008  %select the points above the y=x
    y_x=y_x(index_y)
    Y=[Y,y_x]
end
for j=1:length(Y)
Y{j}=real(Y{j})
end
Y_TM=Y
save("TM.mat","Y_TM")

function f=F(x,y)
f=1-(sqrt(1-1/(y^2-x^2))-(1-1/y^2))^2/(sqrt(1-1/(y^2-x^2))+(1-1/y^2))^2*exp(2j*2*pi*sqrt(y^2-x^2)*1.75)
end
