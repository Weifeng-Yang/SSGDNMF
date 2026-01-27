function [ngmar,R,tlabel]=readfile(i)
  if(i==1) 
   E=load('.\Data\orlraws10P.mat');
   ngmar= double(E.X);
   label= double(E.Y);


    M11=double(ngmar)';
    ngmar=normalize(M11,'range');



   tlabel=double(label);
    R=length(unique(tlabel));

  elseif(i==2)    
   E=load('.\Data\pixraw10P.mat');
   ngmar= double(E.X);
   ngmart=ngmar;
   label= double(E.Y);

    M11=double(ngmar)';
    ngmar=normalize(M11,'range');

   tlabel=double(label);
    R=length(unique(tlabel));




  end



%     dimension(1)=100;
%     dimension(2)=ceil(size(ngmar,2)/4);
%     dimension(end+1)=size(ngmar,2);
%     dimension=[size(ngmar,1),dimension];

end


