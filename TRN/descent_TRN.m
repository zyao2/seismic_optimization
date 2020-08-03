%------------------------------------------------------
% INPUT : integer n dimension                         %
%         real,dimension(n) grad current gradient     %
%         character*4 FLAG communication flag         %
% INPUT/OUTPUT : optim_type optim (data structure)    %
%-----------------------------------------------------%
function optim=descent_TRN(n,grad,optim,FLAG)
 false=0;true=1;
  if(strcmp(optim.CG_phase,'INIT'))
     %-----------------------------------------------------%
     % if optim%CG_phase is INIT, initialize the conjugate %
     % gradient process                                    %
     %-----------------------------------------------------%
     optim.residual(:)=grad(:);
     optim.d(:)=-1.*optim.residual(:);
     optim.Hd(:)=0.;
     optim.descent(:)=0.;     
     optim.qk_CG=0.;
     optim.hessian_term=0.;     
     optim.norm_residual= normL2(n,optim.residual);
     optim.conv_CG=false;    
     optim.cpt_iter_CG=0;
    
     optim.CG_phase='IRUN';
  else                         
     %-----------------------------------------------------%
     % else perform one conjugate gradient iteration       %
     %-----------------------------------------------------%     
     dHd=scalL2(n,optim.d,optim.Hd);        
     if(dHd<0.) 
        %-----------------------------------------------------%
        % if dHd < 0, detection of a negative eigenvalue of   %
        % the Hessian operator, stop the process              %        
        %-----------------------------------------------------%     
        optim.conv_CG=true;
        %write(21,*) 'Negative curvature'
        if(optim.cpt_iter_CG==0)                     
           %-----------------------------------------------------%
           % if this is the first iteration, thenreturn the      %
           % opposite of the gradient as descent direction       %
           % (steepest descent direction)                        %
           %-----------------------------------------------------%     
           optim.descent(:)=optim.d(:);            
           %-----------------------------------------------------%
           % if the debug option is activated, compute the       %
           % quadratic function minimized during the conjugate   %
           % gradient process (check is this value decresae      %
           % throughout the CG iterations )                      %
           %-----------------------------------------------------%     
           if(optim.debug)
              mgrad=zeros(n,1);
              mgrad(:)=-1.*grad(:);              
              optim.norm_residual= normL2(n,optim.residual);
              alpha=(optim.norm_residual.^2)/dHd;         
              optim.qkm1_CG=optim.qk_CG;
              grad_term=scalL2(n,optim.descent,mgrad);          
              optim.hessian_term=optim.hessian_term+(alpha.^2)*dHd;
              optim.qk_CG=-grad_term+0.5*optim.hessian_term;  
              clear mgrad;
           end
        end
     else 
        %-----------------------------------------------------%
        % if dHd > 0, then perform one conjugate gradient     %
        % iteration                                           %
        %-----------------------------------------------------%     
        %Update descent direction
        optim.norm_residual= normL2(n,optim.residual);
        alpha=(optim.norm_residual.^2)/dHd;
        optim.descent_prev(:)=optim.descent(:);
        optim.descent(:)=optim.descent(:)+alpha*optim.d(:);
        optim.residual(:)=optim.residual(:)+alpha*optim.Hd(:);  
        %Update CG direction
        norm_residual_prev=optim.norm_residual;
        optim.norm_residual= normL2(n,optim.residual);
        beta=(optim.norm_residual.^2)/(norm_residual_prev.^2);
        optim.d(:)=-1.*optim.residual(:)+beta*optim.d(:);                
        %Update iteration counter 
        optim.cpt_iter_CG=optim.cpt_iter_CG+1;
        %-----------------------------------------------------%
        % if the debug option is activated, compute the       %
        % quadratic function minimized during the conjugate   %
        % gradient process (check is this value decresae      %
        % throughout the CG iterations )                      %
        %-----------------------------------------------------%        
        if(optim.debug)
           mgrad=zeros(n,1);
           mgrad(:)=-1.*grad(:) ;             
           optim.qkm1_CG=optim.qk_CG;
           grad_term=scalL2(n,optim.descent,mgrad);
           descent_scal_Hd= scalL2(n,optim.descent_prev,optim.Hd);
           optim.hessian_term=optim.hessian_term+(alpha.^2)*dHd+...
                2.*alpha*descent_scal_Hd;
           optim.qk_CG=-grad_term+0.5*optim.hessian_term;
           clear mgrad;
        end
        %-----------------------------------------------------%
        % Check if the Eisenstat stopping critertion is       %
        % satisfied                                           %
        %-----------------------------------------------------%        
        optim.conv_CG=((optim.norm_residual<=(optim.eta*optim.norm_grad))|| ...
             (optim.cpt_iter_CG>=optim.niter_max_CG));       
        %-----------------------------------------------------%
        % Print information on the current CG iteration       %
        %-----------------------------------------------------%
        %call print_info_TRN(n,optim,0e0,FLAG)     
     end
  end
  
end %subroutine descent_TRN