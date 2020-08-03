%-----------------------------------------------------%
% INPUT  : integer :: n (dimension)                   % 
%          real fcost (current cost)                  %
%          real,dimension(n) grad                     %
%          real,dimension(n) grad_preco               %
% INPUT/OUTPUT : real,dimension(n) x                  %
%                optim_typ optim (data structure)     %
%                character*4 FLAG (communication)     %
%-----------------------------------------------------%
function [x,optim,FLAG]= PLBFGS(n,x,fcost,grad,grad_preco,optim,FLAG)

  
  if(strcmp(FLAG,'INIT'))
     %-----------------------------------------------------%
     % if FLAG is INIT, call the dedicated initialization  %
     % subroutine to allocate data structure optim and     %
     % initialize the linesearch process                   %
     %-----------------------------------------------------%
     optim= init_PLBFGS(n,optim.l,x,fcost,grad,grad_preco,optim);
     [x,optim]= std_linesearch(n,x,fcost,grad,optim);
     
     FLAG='GRAD';     
     optim.nfwd_pb=optim.nfwd_pb+1;
  elseif(strcmp(FLAG,'PREC'))
     %-----------------------------------------------------%
     % if FLAG is PREC, we return from a call to the       %
     % user preconditioner, then we have to finish the     %
     % computation of the descent direction                %
     %-----------------------------------------------------%
     optim.descent= descent2_PLBFGS(n,grad,optim.sk,optim.yk,optim.cpt_lbfgs,...
          optim.l,optim);         
     %LBFGS save
     optim= save_LBFGS(n,optim.cpt_lbfgs,optim.l,x,grad,optim);
     optim.cpt_iter=optim.cpt_iter+1;        
     %-----------------------------------------------------%
     % before continuing we test for convergence           %
     %-----------------------------------------------------%
     test_conv=std_test_conv(optim,fcost);
     if(test_conv)
        FLAG='CONV';
        %print info on current iteration
        optim.grad(:)=grad(:);
        
        optim= finalize_PLBFGS(optim);
     else
        FLAG='NSTE';
        %print info on current iteration
        optim.grad(:)=grad(:);
        
     end
  else
     %-----------------------------------------------------%
     % else call the linesearch process                    %  
     %-----------------------------------------------------%
     [x,optim]= std_linesearch(n,x,fcost,grad,optim);   
     if(strcmp(optim.task,'NEW_STEP'))
        %LBFGS update        
        optim= update_LBFGS(n,optim.cpt_lbfgs,optim.l,x,grad,optim);
        %Start the computation of the new descent direction
        optim= descent1_PLBFGS(n,grad,optim.sk,optim.yk,optim.cpt_lbfgs,...
             optim.l,optim);         
        %Set FLAG to PREC for asking user to perform preconditioning 
        FLAG='PREC';                        
     elseif(strcmp(optim.task,'NEW_GRAD'))
        %-----------------------------------------------------%
        % if the linesearch needs a new gradient then ask the %  
        % user to provide it
        %-----------------------------------------------------%
        FLAG='GRAD';         
        optim.nfwd_pb=optim.nfwd_pb+1;
     elseif(strcmp(optim.task,'FAILURE%'))    
        %-----------------------------------------------------%
        % if the linesearch has failed, inform the user       %
        %-----------------------------------------------------%
        FLAG='FAIL';
        %print info on current iteration
        optim.grad(:)=grad(:);
       
        optim= finalize_PLBFGS(optim);
     end
  end
  
end % PLBFGS



