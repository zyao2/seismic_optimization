%*****************************************************%
%*          SEISCOPE OPTIMIZATION TOOLBOX            *%
%*****************************************************%
% This routine is the reverse communication mode      %
% truncated Newton. This algorithm is described in    %
% L.Metivier, R. Brossier, J. Virieux, S.Operto,      %
% Truncated Newton and full waveform inversion, 2013  %
% SIAM Journal on Scientific Computing, Vol. 35,      %
% No. 2, pp. B401–B437,                               %         
%                                                     %
% This routine performs an iterative                  %
% minimization of a function f following the          %
% recurrence                                          %
%                                                     %
% x_{0}=x0                                            %
% x_{k+1}=x_k+\alpha_k d_k  (1)                       %
%                                                     %
% where the descent direction d_k is computed through %
% the resolution of the linear system                 %
%                                                     %
% H_k d_k=- \nabla f_k  (2)                           %
%                                                     %
% with H_k        : Hessian operator at iteration k   %
%                                                     %
%\nabla f_k       : gradient of f in x_k              %
%                                                     %
% and alpha_k is the steplength computed through the  %
% common linesearch algorithm of the TOOLBOX          %
%                                                     %
% The linear system (2) is solved through a matrix    %
% free conjugate gradient algorithm which requires the%
% user to perform multiplication of given vector by   %
% the Hessian operator.                               %                     
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
% - FLAG='GRAD' => the user must compute the cost and %
%                  gradient at current point x in     %
%                  fcost and grad                     %
% - FLAG='HESS' => the user must multiply the vector  %
%                  optim%d by the Hessian operator and%
%                  set the result in optim%Hd         %
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
function [x,optim,FLAG]= TRN(n,x,fcost,grad,optim,FLAG)
% include 'optim_type.h'  
  if(strcmp(FLAG,'INIT'))
     %-----------------------------------------------------%
     % if FLAG is INIT, call the dedicated initialization  %
     % subroutine to allocate data structure optim         %
     %-----------------------------------------------------%
     optim= init_TRN(n,x,fcost,grad,optim);       
     optim.comm='DESC' ;    
     optim.CG_phase='INIT';
     optim.nfwd_pb=optim.nfwd_pb+1;
     optim.conv_CG=0;
     FLAG='NONE';
  end
  if(strcmp(optim.comm,'DESC')) 
     %-----------------------------------------------------%
     % if optim%comm is DES, the optimizer is computing    %
     % a descent direction through the conjugate gradient  %
     %-----------------------------------------------------%
     optim= descent_TRN(n,grad,optim,FLAG);     
     if(optim.conv_CG)       
         %-----------------------------------------------------%
        % if the conjugate gradient has converged go to next  %
        % phase: linesearch in the descent direction          %
        %-----------------------------------------------------%  
        optim.comm='NSTE';
        optim.CG_phase='INIT';        
        FLAG='NONE' ;       
     else
        %-----------------------------------------------------%
        % else perform a new iteration of conjugate gradient  % 
        % and ask the user to compute a Hessian-vector product%
        %-----------------------------------------------------%  
        FLAG='HESS';
        optim.nhess=optim.nhess+1;
     end
  elseif(optim.comm=='NSTE') 
     %-----------------------------------------------------%
     % if optim%comm is NSTE, the optimizer is looking for %
     % a new step in the descent direction                 %
     %-----------------------------------------------------%
     [x,optim]= std_linesearch(n,x,fcost,grad,optim);     
     if(optim.task=='NEW_STEP') %NEW STEP            
        %-----------------------------------------------------%
        % if optim%task is 'NEW_STEP, the linesearch process  %
        % has found the new step                              %
        %-----------------------------------------------------%
        optim.cpt_iter=optim.cpt_iter+1;
        %Save the previous gradient norm 
        optim.norm_grad_m1=optim.norm_grad;
        %Compute the new gradient norm 
        optim.norm_grad= normL2(n,grad);        
        %Print info on current nonlinear iteration
       % call print_info_TRN(n,optim,fcost,FLAG);        
        %Test for convergence
        test_conv= std_test_conv(optim,fcost);
        if(test_conv) 
           FLAG='CONV';          
          % call print_info_TRN(n,optim,fcost,FLAG) ;                 
           optim= finalize_TRN(optim);
        else
           %Flags for the computation of the new descent direction
           FLAG='NSTE';           
           optim.comm='DESC';       
           %Update forcing term optim%eta following the Eisenstat and Walker formula
           optim= forcing_term_TRN(n,grad,optim) ;       
        end
     elseif(optim.task=='NEW_GRAD')  %STILL SEARCHING THE STEP 
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
        %call print_info_TRN(n,optim,fcost,FLAG);

        optim=finalize_TRN(optim);
     end
  end
end % TRN