clearvars 
clc

%% Parameter.
%   dimension : The ranks of matrix decomposition
%   index     : The dataset to be used, when index=1, use orlraws10P dataset
%               when index=2, use pixraw10P dataset. 
%   r         : Step factor
%   maxiteropt: Maximum iteration alloted to the method
%   trigger   : Whether to enable the indicator array of each method, where
%               when 1∈trigger, enable the 𝓁21-NMF method
%               when 2∈trigger, enable the 𝓁21-GSNMF method
%               when 3∈trigger, enable the GDNMF method
%               when 4∈trigger, enable the 𝓁1-SDNMF method
%               when 5∈trigger, enable the DHGLpSNMF method
%               when 6∈trigger, enable the 𝓁1-SGDNMF method 
%               when 7∈trigger, enable the 𝓁21-SSDNMF method
%               when 8∈trigger, enable the 𝓁21-SSGDNMF method
%   lamda     ：The sparse regularization parameters
%   alphat    : The graph regularization parameter
%   stopindex : The indicator of the stop condition.  
%               To set the specific termination condition, see the 'stopcheck' function for details.  
%               The default termination condition is: ϵ<1e-6 or maxiteropt>8000 
%% Display
%   nonzero       ：The number of non-zero elements in each component.
%   nonzero rows  ：The number of non-zero rows in each component.
%   Rel           : The difference in the variable value between two iterations.


%% Parameter settings and Select dataset
rng('shuffle')
warning('off');
index=1;
maxiteropt=8000;
r=1.01;
rho=10^-8;
trigger=[1,2,3,4,5,6,7,8];
lamda=[1,1,1];
alphat=1;
[ngmar,N,y]=readfile(index);
dimension=[size(ngmar,1),100,ceil(size(ngmar,2)/4),size(ngmar,2)];
stopindex=4;




    


num=length(dimension)-1;
for j=1:10
%% Init
for i=1:num
    var{i}=rand(dimension(i),dimension(i+1));
end

%% Solving
for i=1:length(trigger)       
[datas{i},vars{i}]=ALGOchoose(var,ngmar,maxiteropt,trigger(i),stopindex,r,rho,alphat,lamda);
end
datas{length(trigger)+1}=vars;
datass{j}=datas;

for i=1:length(trigger)
    vars1=datas{length(trigger)+1};
    vartemp=vars1{i};
    [acc(j,i),rdx(j,i),NMIs(j,i)]=clustermeans(vartemp{end}',N,y);
end


end






 
 
%% Display results
accmean=mean(acc)
rdxmean=mean(rdx)
nmimean=mean(NMIs)
plt0=barplot(trigger,acc,rdx,NMIs); plt0=barplot(trigger,acc,rdx,NMIs);

 
 
 



