%-----------------------------------------------------%
% INPUT : integer n dimension                         %
%         integer cpt_lbfgs counts the number of      % 
%                           already stored pairs      %
%         integer l maximum number of stored pairs    %
%         real,dimension(n) x current point           %
%         real,dimension(n) x current gradient        %
% OUTPUT : sk,yk updated vector                       %
%-----------------------------------------------------%
function optim=save_LBFGS(n,cpt_lbfgs,l,x,grad,optim)

  
  if(cpt_lbfgs<=l)
     %---------------------------------------------------%
     % if the number of stored pairs does not exceed the %
     % maximum value, then save x and grad               %
     %---------------------------------------------------%
     optim.sk(:,cpt_lbfgs)=x(:);
     optim.yk(:,cpt_lbfgs)=grad(:);   
  else
     %---------------------------------------------------%
     % otherwise, erase the oldest pair and save the     %
     % new one (shift)                                   %
     %---------------------------------------------------%
     for i=1:l-1
        optim.sk(:,i)=optim.sk(:,i+1);
        optim.yk(:,i)=optim.yk(:,i+1);                
     end
     optim.sk(:,l)=x(:);
     optim.yk(:,l)=grad(:);
  end
  
end % save_LBFGS