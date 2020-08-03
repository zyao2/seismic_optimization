% This routine performs an iterative                  %
% minimization of a function f following the          %
% recurrence                                          %
%                                                     %
% x_{0}=x0                                            %
% x_{k+1}=x_k+\alpha_k d_k                            %
%                                                     %
% where the descent direction d_k is                  %
%                                                     %
% d_k=-Q_k \nabla f_k + \beta_k d_{k-1}               %
%                                                     %
% where Q_k       : preconditioner at iteration k     %
%      \nabla f_k : gradient of f in x_k              %
%      \beta_k    : scalar computed through the       %
%                   Dai and Yuan formula              %
%                                                     %
% and alpha_k is the steplength computed through the  %
% common linesearch algorithm of the TOOLBOX          %
%                                                     %
% The first call to the algorithm must be done with   %
% FLAG='INIT'. For this first call, the initial point %
% x0 is given through the variable x. The input       %
% variables fcost and grad, grad_preco must correspond%
% respectively to the misfit, gradient and            % 
% preconditioned gradient at x0.                      %
%                                                     %
% The reverse communication with the user is          %
% performed through the variable FLAG. This           %
% variable indicates to the user on return what action% 
% he has to do, or the state of the algorithm.        %
% Possible values are                                 %
% - FLAG='GRAD' => the user must compute the cost,    %
%                  gradient and preconditioned        %
%                  gradient at current point x        % 
% - FLAG='CONV' => a minimizer has been found         %
% - FLAG='NSTE' => a new step is performed            %
% - FLAG='FAIL' => the linesearch has failed          %
%-----------------------------------------------------%
% INPUT  : integer :: n (dimension)                   % 
%          real fcost (current cost)                  %
%          real,dimension(n) grad                     %
%          real,dimension(n) grad_preco               %
% INPUT/OUTPUT : real,dimension(n) x                  %
%                optim_typ optim (data structure)     %
%                character*4 FLAG (communication)     %
%-----------------------------------------------------%
function [x,optim,FLAG]= PNLCG(n,x,fcost,grad,grad_preco,optim,FLAG)

  
  if(strcmp(FLAG,'INIT'))
     %-----------------------------------------------------%
     % if FLAG is INIT, call the dedicated initialization  %
     % subroutine to allocate data structure optim and     %
     % initialize the linesearch process                   %
     %-----------------------------------------------------%
     optim= init_PNLCG(n,x,fcost,grad,grad_preco,optim);
     [x,optim]=std_linesearch(n,x,fcost,grad,optim);
    
     FLAG='GRAD' ;    
     optim.nfwd_pb=optim.nfwd_pb+1;
     %Store current gradient before the user compute the new one
     optim.grad_prev(:)=grad(:);      
  else
     %-----------------------------------------------------%
     % else call the linesearch process                    %  
     %-----------------------------------------------------%
     [x,optim]=std_linesearch(n,x,fcost,grad,optim);     
     if(strcmp(optim.task,'NEW_STEP'))
        optim.cpt_iter=optim.cpt_iter+1;
        %-----------------------------------------------------%
        % test for convergence                                %
        %-----------------------------------------------------%
        test_conv= std_test_conv(optim,fcost);
        if(test_conv)
           FLAG='CONV';
           %print info on current iteration        
           optim.grad(:)=grad(:);
          
           optim= finalize_PNLCG(optim);
        else
           FLAG='NSTE';
           %-----------------------------------------------------%
           % if a NEW_STEP is taken, compute a new descent       %
           % direction using current descent, gradient and       %
           % preconditioned gradient                             %
           %-----------------------------------------------------%
           optim= descent_PNLCG(n,grad,grad_preco,optim);                    
           %print info on current iteration        
           optim.grad(:)=grad(:);
           
        end
     elseif(strcmp(optim.task,'NEW_GRAD'))
        %-----------------------------------------------------%
        % if the linesearch needs a new gradient then ask the %  
        % user to provide it
        %-----------------------------------------------------%
        FLAG='GRAD';         
        optim.nfwd_pb=optim.nfwd_pb+1;
        %Store current gradient before the user compute the new one
        optim.grad_prev(:)=grad(:);        
     elseif(strcmp(optim.task,'FAILURE%'))      
        %-----------------------------------------------------%
        % if the linesearch has failed, inform the user       %
        %-----------------------------------------------------%
        FLAG='FAIL';
        %print info on current iteration
        optim.grad(:)=grad(:);
        
        optim= finalize_PNLCG(optim);
     end
  end
  
end %subroutine PNLCG