%-----------------------------------------------------%
% INPUT : integer n dimension                         %
%         integer cpt_lbfgs counts the number of      % 
%                           already stored pairs      %
%         integer l maximum number of stored pairs    %
%         real,dimension(n) x current point           %
%         real,dimension(n) x current gradient        %
% OUTPUT : sk,yk updated vectors                      %
%-----------------------------------------------------%
function optim=update_LBFGS(n,cpt_lbfgs,l,x,grad,optim)


  if(cpt_lbfgs<=l) 
     %---------------------------------------------------%
     % if the number of stored pairs does not exceed the %
     % maximum value, then compute a new pair sk yk and  %
     % update the counter cpt_lbfgs                      %
     %---------------------------------------------------%
     optim.sk(:,cpt_lbfgs)=x(:)-optim.sk(:,cpt_lbfgs);
     optim.yk(:,cpt_lbfgs)=grad(:)-optim.yk(:,cpt_lbfgs);
     optim.cpt_lbfgs=optim.cpt_lbfgs+1;
  else
     %---------------------------------------------------%
     % otherwise, simply update the lth pair             %
     %---------------------------------------------------%
     optim.sk(:,l)=x(:)-optim.sk(:,l);
     optim.yk(:,l)=grad(:)-optim.yk(:,l);
  end
  
end % update_LBFGS
