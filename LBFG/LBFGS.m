%-----------------------------------------------------%
% INPUT  : integer :: n (dimension)                   % 
%          real fcost (current cost)                  %
%          real,dimension(n) grad                     %
% INPUT/OUTPUT : real,dimension(n) x                  %
%                optim_type optim (data structure)    %
%                character*4 FLAG (communication)     %
%-----------------------------------------------------%
function [x,optim,FLAG]=LBFGS(n,x,fcost,grad,optim,FLAG)

   if(FLAG=='INIT')
     %-----------------------------------------------------%
     % if FLAG is INIT, call the dedicated initialization  %
     % subroutine to allocate data structure optim and     %
     % initialize the linesearch process                   %
     %-----------------------------------------------------%
     optim= init_LBFGS(n,optim.l,x,fcost,grad,optim);
     [x,optim]= std_linesearch(n,x,fcost,grad,optim);
     FLAG='GRAD';     
     optim.nfwd_pb=optim.nfwd_pb+1;
  else
     %-----------------------------------------------------%
     % else call the linesearch process                    %  
     %-----------------------------------------------------%
     [x,optim]= std_linesearch(n,x,fcost,grad,optim); 
     if(strcmp(optim.task,'NEW_STEP'))
        %-----------------------------------------------------%
        % test for convergence                                %
        %-----------------------------------------------------%
        optim.cpt_iter=optim.cpt_iter+1;                   
        test_conv= std_test_conv(optim,fcost);
        if(test_conv) 
           FLAG='CONV';           
           %print info on current iteration
           optim.grad(:)=grad(:);
        else
           FLAG='NSTE';
           %-----------------------------------------------------%
           % if a NEW_STEP is taken, compute a new descent       %
           % direction using current gradient and l-BFGS         %
           % approximation of the inverse Hessian                %
           % preconditioned gradient                             %
           %-----------------------------------------------------%
           %LBFGS update        
           %[optim.sk,optim.yk]= update_LBFGS(n,optim.cpt_lbfgs,optim.l,x,grad,optim.sk,optim.yk);
           optim= update_LBFGS(n,optim.cpt_lbfgs,optim.l,x,grad,optim);
           %Computation of the new descent direction
           optim.descent= descent_LBFGS(n,grad,optim.sk,optim.yk,optim.cpt_lbfgs,...
                optim.l);         
           %LBFGS store
           optim= save_LBFGS(n,optim.cpt_lbfgs,optim.l,x,grad,optim);           
           %print info on current iteration
           optim.grad(:)=grad(:);
           
        end
     elseif(optim.task=='NEW_GRAD')
        %-----------------------------------------------------%
        % if the linesearch needs a new gradient then ask the %  
        % user to provide it
        %-----------------------------------------------------%
        FLAG='GRAD';         
        optim.nfwd_pb=optim.nfwd_pb+1;        
     elseif(optim.task=='FAILURE%')       
        %-----------------------------------------------------%
        % if the linesearch has failed, inform the user       %
        %-----------------------------------------------------%
        FLAG='FAIL';
        %print info on current iteration
        optim.grad(:)=grad(:);
        
        optim= finalize_LBFGS(optim);
     end
  end
  
end %LBFGS
