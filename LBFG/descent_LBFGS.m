%-----------------------------------------------------%
% INPUT : integer n dimension                         %
%         real,dimension(n) grad current gradient     %
%         real,dimension(n,l) sk yk pairs of vectors  %
%                             for l-BFGS approximation%
%         integer cpt_lbfgs current number of stored  %
%                           pairs                     %
%         integer l maximum number of stored pairs    %
% OUTPUT : real,dimension(n) descent                  %
%-----------------------------------------------------%
function descent=descent_LBFGS(n,grad,sk,yk,cpt_lbfgs,l)

  borne_i=cpt_lbfgs-1;
  descent=zeros(n,1);
  %------------------------------------------%
  % SAFEGUARD                                %
  %------------------------------------------%
  norml2_sk= normL2(n,sk(:,borne_i));
  norml2_yk= normL2(n,yk(:,borne_i));
  if( (norml2_sk==0.)||(norml2_yk==0.))
     descent(:)=-1.*grad(:);     
  else
     %------------------------------------------%
     % First phase of the recursion loop        %
     %------------------------------------------%
     alpha=zeros(cpt_lbfgs,1);
     rho=zeros(cpt_lbfgs,1);
     q=zeros(n,1);
     q(:)=grad(:);
     for i=1:borne_i
        rho(borne_i-i+1)= scalL2(n,yk(:,borne_i-i+1),sk(:,borne_i-i+1));
        rho(borne_i-i+1)=1./rho(borne_i-i+1);
        alpha(borne_i-i+1)= scalL2(n,sk(:,borne_i-i+1),q(:));
        alpha(borne_i-i+1)=rho(borne_i-i+1)*alpha(borne_i-i+1);
        q(:)=q(:)-alpha(borne_i-i+1)*yk(:,borne_i-i+1);
     end
     gamma_num= scalL2(n,sk(:,borne_i),yk(:,borne_i));
     gamma_den= normL2(n,yk(:,borne_i));
     %------------------------------------------%
     % Scaling by gamma                         %
     %------------------------------------------%
     gamma=gamma_num/(gamma_den.^2);     
     descent(:)=gamma*q(:);
     %------------------------------------------%
     % Second phase of the recursion loop       %
     %------------------------------------------%
     for i=1:borne_i
        beta= scalL2(n,yk(:,i),descent(:));
        beta=rho(i)*beta;
        descent(:)=descent(:)+(alpha(i)-beta)*sk(:,i);
     end
     descent(:)=-1.*descent(:);
     
  end
  
end% descent_LBFGS