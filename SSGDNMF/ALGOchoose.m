%% Input.
% var         : Initial matrix
% ngmar       : Decomposed matrix
% The remaining parameters are explained the same as the 'main_Run_me' function

%% Output.
% vars:       : Decomposition matrix resulting from the final iterative result
% loss:       : Array of loss functions generated during iteration
% tr:         : Runtime array during iteration


%% An Inertial Block Proximal Linearized Method with Adaptive Momentum for Nonconvex and Nonsmooth Optimization

function [data,varss]=ALGOchoose(var,ngmar,maxiteropt,flag,stopindex,r,rho,alphat,lamda)

if(flag==1)
[vars,loss,tr]=L21NMF(var,ngmar,maxiteropt,stopindex,r,2);
varss=vars;
lossdata=loss;
trdata=tr;


elseif(flag==2)
[vars,loss,tr]=L21GSNMF(var,ngmar,maxiteropt,stopindex,r,2);
varss=vars;
lossdata=loss;
trdata=tr;

elseif(flag==3) 
[vars,loss,tr]=GDNMF(var,ngmar,maxiteropt,stopindex,r,0.9999,rho,0);  
varss=vars;
lossdata=loss;
trdata=tr;


elseif(flag==4) 
[vars,loss,tr]=SDNMF(var,ngmar,maxiteropt,stopindex,r,0.9999,rho,1);
varss=vars;
lossdata=loss;
trdata=tr;





elseif(flag==5) 
[vars,loss,tr]=DHGLpSNMF(var,ngmar,maxiteropt,stopindex,r,0.005,0.8);
varss=vars;
lossdata=loss;
trdata=tr;

elseif(flag==6) 
[vars,loss,tr]=L1SGDNMF(var,ngmar,maxiteropt,stopindex,r,0.9999,rho,1,1);  
varss=vars;
lossdata=loss;
trdata=tr;

elseif(flag==7) 
[vars,loss,tr]=L21SSGDNMF(var,ngmar,maxiteropt,stopindex,r,0.9999,rho,0,lamda);  
varss=vars;
lossdata=loss;
trdata=tr;




elseif(flag==8) 
[vars,loss,tr]=L21SSGDNMF(var,ngmar,maxiteropt,stopindex,r,0.9999,rho,alphat,lamda);  
varss=vars;
lossdata=loss;
trdata=tr;






end











data{1}=lossdata;
data{2}=trdata;


end



