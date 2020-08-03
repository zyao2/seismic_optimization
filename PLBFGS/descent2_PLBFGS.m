%-----------------------------------------------------%
%SECOND LOOP: descent2_PLBFGS                         %
%------------------------------------------------------
% INPUT : integer n dimension                         %
%         real,dimension(n) grad current gradient     %
%         real,dimension(n,l) sk yk pairs of vectors  %
%                             for l-BFGS approximation%
%         integer cpt_lbfgs current number of stored  %
%                           pairs                     %
%         integer l maximum number of stored pairs    %
% OUTPUT : real,dimension(n) descent                  %
%-----------------------------------------------------%
function descent=descent2_PLBFGS(n,grad,sk,yk,cpt_lbfgs,l,optim)
descent=zeros(n,1);
  borne_i=cpt_lbfgs-1;  
  gamma_num=scalL2(n,sk(:,borne_i),yk(:,borne_i));
  gamma_den= normL2(n,yk(:,borne_i));
  gamma=gamma_num/(gamma_den.^2) ;    
  descent(:)=gamma*optim.q_plb(:);
  for i=1:borne_i
     beta= scalL2(n,yk(:,i),descent(:));
     beta=optim.rho_plb(i)*beta;
     descent(:)=descent(:)+(optim.alpha_plb(i)-beta)*sk(:,i);
  end
  descent(:)=-1.*descent(:);
  clear optim.alpha_plb;
  clear optim.rho_plb;
  
end % descent2_PLBFGS
