% Possible values are                                 %
% - FLAG='GRAD' => the user must compute the cost and %
%                  (preconditioned) gradient at       %
%                  current point x                    %
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
function [x,optim,FLAG]= PSTD(n,x,fcost,grad,grad_preco,optim,FLAG)


  if(strcmp(FLAG,'INIT'))    
     %-----------------------------------------------------%
     % if FLAG is INIT, call the dedicated initialization  %
     % subroutine to allocate data structure optim and     %
     % initialize the linesearch process                   %
     %-----------------------------------------------------%
     optim= init_PSTD(n,x,fcost,grad,grad_preco,optim);
     [x,optim]= std_linesearch(n,x,fcost,grad,optim);
     
     FLAG='GRAD' ;    
     optim.nfwd_pb=optim.nfwd_pb+1;
  else
     %-----------------------------------------------------%
     % else call the linesearch process                    %  
     %-----------------------------------------------------%
     [x,optim]=std_linesearch(n,x,fcost,grad,optim);
     if(strcmp(optim.task,'NEW_STEP')) %NEW STEP
        %-----------------------------------------------------%
        % test for convergence                                %
        %-----------------------------------------------------%        
        optim.cpt_iter=optim.cpt_iter+1 ;       
        test_conv= std_test_conv(optim,fcost);
        if(test_conv) 
           FLAG='CONV';

           optim= finalize_PSTD(optim);
        else
           FLAG='NSTE';
           %-----------------------------------------------------%
           % if a NEW_STEP is taken, compute a new descent       %
           % direction using current descent, gradient and       %
           % preconditioned gradient                             %          
           %-----------------------------------------------------%
           optim.grad(:)=grad(:);
           optim.descent(:)=-1.*grad_preco(:);                
 
        end
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
        
        optim= finalize_PSTD(optim);
     end
  end
  
end % PSTD