%%  All parameters of this function are explained the same as 'main_Run_me' and 'ALGOchoose' functions
function [var,loss,timerun]=DHGLpSNMF(var,ngmar,maxiteropt,stopindex,r,beta,p)
%% initialization algorithm
loss=[];
timerun=[0];
num=length(var);


Lapk=LLaplace(ngmar);
LK=zeros(1,num);
L=ones(1,num);
varK=var;
returnloss=norm(ngmar,"fro");


for i=1:num-1
    alpha(i)=0;
end
alpha(num)=0.002;

loss(1)=compute(var,ngmar,beta,p)+alpha(num)/2*trace(var{num}*Lapk{end}*var{num}');


t1=clock;


for i=1:maxiteropt
%% update parameters
fprintf("%d\n",i);
varK=var;
for j=1:num
    if(j==1)
    [V,L(j)]=grad(var,ngmar,Lapk,j,num,r,alpha,beta,p);
%     var{j}=V/norm(V,'fro')^2;
    var{j}=V;
    else
    [V,L(j)]=grad(var,ngmar,Lapk,j,num,r,alpha,beta,p);
    var{j}=PROXn1(V);    
    end
end
loss(i+1)=compute(var,ngmar,beta,p)+alpha(num)/2*trace(var{num}*Lapk{end}*var{num}');




%% Check if termination condition is met
fprintf("DHGLpSNMF\n");
check1=0;
check2=0;
for j=1:num
    fprintf("nonzero:%d\n",nnz(var{j}));
    check1=check1+norm(var{j}-varK{j},'fro');
    check2=check2+norm(varK{j},'fro');
end
t2=clock;
timerun(i+1)=etime(t2,t1);
Res=check1/check2;
fprintf("Rel：%d\n",Res);
stop=stopcheck(Res,timerun,stopindex);
if(stop==1)
    fprintf("Number of terminations：%d\n",i);
    pause(2);
    break;
end





end

end











function x=PROXn1(x)
x=double(x);
x(x<0)=0;
end


function loss=compute(var,ngmar,beta,p)
    nga=var{1};
    for i=2:length(var)
        nga=nga*var{i};
    end
    temp=abs(var{end}).^p;
    loss=norm(ngmar-nga,'fro')+beta*sum(temp(:));
end


%% Gradient calculation function for DNMF
function [U,L]=grad(var,ngmar,Lapk,n,num,r,alpha,beta,p)
    if(n>1)
    Z=var{1};
    end
    if(n<num)
    ZT=var{n+1};
    end

   
    for i=2:n-1
        Z=Z*var{i};
    end

    
    for i=n+2:num
        ZT=ZT*var{i};
    end 
%     toc


    if(n==1)
    ck=norm(ZT*ZT',2);
    ck=checkck(ck);
    L=ck;
    mar=var{n}*ZT;
    U=var{n}-1/(r*L)*(mar-ngmar)*ZT'; 
    elseif(n<num)
    ck=norm(Z'*Z,2)*norm(ZT*ZT',2);
    ck=checkck(ck);
    L=ck;
    mar=Z*var{n}*ZT;
    U=var{n}-1/(r*L)*Z'*(mar-ngmar)*(ZT)';  
    elseif(n==num)
    ck=norm(Z'*Z,2);
    L=ck+alpha(end)*norm(Lapk{end},2)+(beta*p)^(1/(p-1));
    mar=Z*var{n};
    U=var{n}-1/(r*L)*(Z'*(mar-ngmar)+alpha(n)*var{n}*Lapk{end}+beta*p*var{end}.^(p-1));   
    end
end





function ck=checkck(ck)
    if(ck==0) 
        a=0;
       while(a==0)
           a=rand(1);
       end
       ck=a;
     end
end

