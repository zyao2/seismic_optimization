%-----------------------------------------------------%
%FIRST LOOP: descent1_PLBFGS                          %
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
function optim= descent1_PLBFGS(n,grad,sk,yk,cpt_lbfgs,l,optim)

  borne_i=cpt_lbfgs-1;

  optim.alpha_plb=zeros(cpt_lbfgs,1);
  optim.rho_plb=zeros(cpt_lbfgs,1);
  optim.q_plb(:)=grad(:);
  for i=1:borne_i
     optim.rho_plb(borne_i-i+1)= scalL2(n,yk(:,borne_i-i+1),sk(:,borne_i-i+1));
     optim.rho_plb(borne_i-i+1)=1./optim.rho_plb(borne_i-i+1);
     optim.alpha_plb(borne_i-i+1)= scalL2(n,sk(:,borne_i-i+1),optim.q_plb(:));
     optim.alpha_plb(borne_i-i+1)=optim.rho_plb(borne_i-i+1)*optim.alpha_plb(borne_i-i+1);
     optim.q_plb(:)=optim.q_plb(:)-optim.alpha_plb(borne_i-i+1)*yk(:,borne_i-i+1);
  end
  
end % descent1_PLBFGS

