% This routine performs an iterative                  %
% minimization of a function f following the          %
% recurrence                                          %
%                                                     %
% x_{0}=x0                                            %
% x_{k+1}=x_k+\alpha_k d_k  (1)                       %
%                                                     %
% where the descent direction d_k is computed through %
% the resolution of the preconditioned linear system  %
%                                                     %
% P_k H_k d_k=- P_k\nabla f_k  (2)                    %
%                                                     %
% with H_k        : Hessian operator at iteration k   %
%                                                     %
% with P_k        : Preconditioner                    %
%                                                     %
%\nabla f_k       : gradient of f in x_k              %
%                                                     %
% and alpha_k is the steplength computed through the  %
% common linesearch algorithm of the TOOLBOX          %
%                                                     %
% The linear system (2) is solved through a matrix    %
% free conjugate gradient algorithm which requires the%
% user to perform multiplication of given vector by   %
% the Hessian operator and the preconditioner.        %                     
%                                                     %
% Because of these tow nested algorithms, the reverse %
% communication strategy requires additional          %
% communicators within the code to clearly track which%
% stage the optimizer has reached                     %
%                                                     %
% The first call to the algorithm must be done with   %
% FLAG='INIT'. For this first call, the initial point %
% x0 is given through the variable x, and the input   %
% variable fcost and grad must correspond respectively%
% to the misfit and gradient at x0.                   %
%                                                     %
% The reverse communication with the user is          %
% performed through the variable FLAG. This           %
% variable indicates to the user on return what action% 
% he has to do, or the state of the algorithm.        %
% Possible values are                                 %
% - FLAG='GRAD' => the user must compute the cost,the %
%                  gradient and the preconditioned    %
%                  gradient at current point x in     %
%                  fcost, grad, grad_preco            %
% - FLAG='HESS' => the user must multiply the vector  %
%                  optim%d by the Hessian operator and%
%                  set the result in optim%Hd         %
% - FLAG='PREC' => the user must multiply the vector  %
%                  optim%residual by the              %
%                  preconditioner and set the result  %
%                  in optim%residual_preco            %
% - FLAG='CONV' => a minimizer has been found         %
% - FLAG='NSTE' => a new step is performed            %
% - FLAG='FAIL' => the linesearch has failed          %
%-----------------------------------------------------%
% INPUT  : integer :: n (dimension)                   % 
%          real fcost (current cost)                  %
%          real,dimension(n) grad                     %
% INPUT/OUTPUT : real,dimension(n) x                  %
%                optim_typ optim (data structure)     %
%                character*4 FLAG (communication)     %
%-----------------------------------------------------%
function [x,optim,FLAG]= PTRN(n,x,fcost,grad,grad_preco,optim,FLAG)
false=0;
  if(strcmp(FLAG,'INIT'))
     %-----------------------------------------------------%
     % if FLAG is INIT, call the dedicated initialization  %
     % subroutine to allocate data structure optim         %
     %-----------------------------------------------------%
     optim= init_PTRN(n,x,fcost,grad,optim);     
     % print_info_PTRN(optim,fcost,FLAG)     
      optim.comm='DES1';     
     optim.CG_phase='INIT';
     optim.nfwd_pb=optim.nfwd_pb+1;
     optim.conv_CG=false;
     FLAG='NONE';
  end
  if(strcmp(optim.comm,'DES1')) 
     %-----------------------------------------------------%
     % if optim%comm is DES1, the optimizer starts the     %
     % computation of a descent direction through the      %
     % conjugate gradient                                  %
     %-----------------------------------------------------%
     if(strcmp(optim.CG_phase,'INIT') )
        %-----------------------------------------------------%
        % if optim%CG_phase is INIT, initialization of the    %
        % conjugate gradient process
        %-----------------------------------------------------%
        optim= descent_PTRN0(n,grad,grad_preco,optim,FLAG);  
        optim.CG_phase='IRUN';
        optim.comm='DES1';
        FLAG='HESS';
        optim.nhess=optim.nhess+1;
      elseif(optim.CG_phase=='IRUN')
        %-----------------------------------------------------%
        % if optim%CG_phase is IRUN, iterate the conjugate    %
        % gradient process (first part)                       %
        %-----------------------------------------------------%
        optim= descent_PTRN1(n,grad,optim,FLAG);    
        if(optim.conv_CG)
           %-----------------------------------------------------%
           % if the conjugate gradient has already converged     % 
           % (detection of a negative curvature), go to next     %
           % phase: linesearch in the descent direction          %
           %-----------------------------------------------------%
           optim.comm='NSTE';
           optim.CG_phase='INIT';        
           FLAG='NONE';                   
        else
           %-----------------------------------------------------%
           % if the conjugate gradient has not converged         % 
           % prematurly, then ask the user to apply the          %
           % preconditioner                                      %
           %-----------------------------------------------------%
           FLAG='PREC';
           optim.comm='DES2';
        end
     end
  elseif(optim.comm=='DES2')
     %-----------------------------------------------------%
     % if optim%comm is DES2, the optimizer finish the     %
     % current iteration of the conjugate gradient process %
     %-----------------------------------------------------%
     optim=descent_PTRN2(n,grad,optim,FLAG);            
     if(optim.conv_CG)
        %-----------------------------------------------------%
        % if the conjugate gradient has converged go to next  %
        % phase: linesearch in the descent direction          %
        %-----------------------------------------------------%        
        optim.comm='NSTE';
        optim.CG_phase='INIT' ;       
        FLAG='NONE';                
     else
        %-----------------------------------------------------%
        % else start a new iteration of conjugate gradient and% 
        % ask the user to compute a Hessian-vector product    %
        %-----------------------------------------------------%        
        optim.comm='DES1';
        FLAG='HESS';
        optim.nhess=optim.nhess+1;
     end
  elseif(optim.comm=='NSTE')         
     %-----------------------------------------------------%
     % if optim%comm is NSTE, a descent direction has been %
     % computed, and a linesearch must be performed in this%
     % direction                                           %
     %-----------------------------------------------------%
      [x,optim]= std_linesearch(n,x,fcost,grad,optim);     
     if(optim.task=='NEW_STEP')  %NEW STEP                  
        %-----------------------------------------------------%
        % if optim%task is 'NEW_STEP, the linesearch process  %
        % has found the new step                              %
        %-----------------------------------------------------%
        optim.cpt_iter=optim.cpt_iter+1;
        %Save the previous gradient norm 
        optim.norm_grad_m1=optim.norm_grad;
        %Computation of the new gradient norm 
        optim.norm_grad= normL2(n,grad);
        %Print infor on current nonlinear iteration
        % print_info_PTRN(optim,fcost,FLAG);        
        %Test for convergence
        test_conv=std_test_conv(optim,fcost);
        if(test_conv) 
           FLAG='CONV';
           % print_info_PTRN(optim,fcost,FLAG)        
           
           optim= finalize_PTRN(optim);
        else
           FLAG='NSTE';
           %Flag for the computation of the new descent direction
           optim.comm='DES1';       
           %Update forcing term optim%eta following the Eisenstat and Walker formula
           optim= forcing_term_TRN(n,grad,optim) ; 
        end
     elseif(optim.task=='NEW_GRAD')
        %-----------------------------------------------------%
        % if optim%task is 'NEW_GRAD, the linesearch process  %
        % is continuing, the gradient at the current point is %
        % required                                            % 
        %-----------------------------------------------------%
        FLAG='GRAD';         
        optim.nfwd_pb=optim.nfwd_pb+1;
     elseif(optim.task=='FAILURE')     
        %-----------------------------------------------------%
        % if optim%task is 'FAILURE, the linesearch process   %
        % has failed, the iterations are stopped              %
        %-----------------------------------------------------%
        FLAG='FAIL';
        %print_info_PTRN(optim,fcost,FLAG)
        optim= finalize_PTRN(optim);
     end
  end
  
end % PTRN



