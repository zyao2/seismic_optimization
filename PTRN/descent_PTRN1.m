%-----------------------------------------------------%
% descent_PTRN1
%-----------------------------------------------------%
% INPUT : integer n dimension                         %
%         real,dimension(n) grad current gradient     %
%         character*4 FLAG communication flag         %
% INPUT/OUTPUT : optim_type optim (data structure)    %
%-----------------------------------------------------%
function optim= descent_PTRN1(n,grad,optim,FLAG)
 true=1;
  %-----------------------------------------------------%
  % start one conjugate gradient iteration       %
  %-----------------------------------------------------%    
  optim.dHd= scalL2(n,optim.d,optim.Hd);        
  if(optim.dHd<0.)
     %-----------------------------------------------------%
     % if dHd < 0, detection of a negative eigenvalue of   %
     % the Hessian operator, stop the process              %        
     %-----------------------------------------------------%  
     optim.conv_CG=true;
     if(optim.cpt_iter_CG==0)
        %-----------------------------------------------------%
        % if this is the first iteration, thenreturn the      %
        % opposite of the preconditioned gradient as descent  %
        % direction: preconditioned steepest descent direction%
        %-----------------------------------------------------%   
        optim.descent(:)=optim.d(:) ;
        %-----------------------------------------------------%
        % if the debug option is activated, compute the       %
        % quadratic function minimized during the conjugate   %
        % gradient process (check is this value decresae      %
        % throughout the CG iterations )                      %
        %-----------------------------------------------------%
        if(optim.debug)
           optim.res_scal_respreco= scalL2(n,optim.residual,optim.residual_preco);
           optim.alpha_CG=optim.res_scal_respreco/optim.dHd;         
           optim.qkm1_CG=optim.qk_CG;
           mgrad=zeros(n,1);
           mgrad(:)=-1.*grad(:) ;          
           grad_term= scalL2(n,optim.descent,mgrad);           
           optim.hessian_term=optim.hessian_term+(optim.alpha_CG.^2)*optim.dHd;
           optim.qk_CG=-grad_term+0.5*optim.hessian_term; 
           clear mgrad;
        end
     end
  else         
     %-----------------------------------------------------%
     % if dHd > 0, then start one conjugate gradient       %
     % iteration                                           %
     %-----------------------------------------------------%   
     %Update descent direction
     optim.res_scal_respreco= scalL2(n,optim.residual,optim.residual_preco);
     optim.alpha_CG=optim.res_scal_respreco/optim.dHd;
     optim.descent_prev(:)=optim.descent(:);
     optim.descent(:)=optim.descent(:)+optim.alpha_CG*optim.d(:);
     optim.residual(:)=optim.residual(:)+optim.alpha_CG*optim.Hd(:);       
     %STOP HERE and wait for preconditioning
  end 
  
end %descent_PTRN1


