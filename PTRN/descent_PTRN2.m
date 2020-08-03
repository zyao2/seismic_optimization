%-----------------------------------------------------%
% descent_PTRN2
%-----------------------------------------------------%
% INPUT : integer n dimension                         %
%         real,dimension(n) grad current gradient     %
%         character*4 FLAG communication flag         %
% INPUT/OUTPUT : optim_type optim (data structure)    %
%-----------------------------------------------------%
function optim= descent_PTRN2(n,grad,optim,FLAG)

  %-----------------------------------------------------%
  % continue the current conjugate gradient iteration   %
  %-----------------------------------------------------% 
  %Update CG direction
  res_scal_respreco_prev=optim.res_scal_respreco;
  optim.res_scal_respreco= scalL2(n,optim.residual,optim.residual_preco);
  beta=(optim.res_scal_respreco)/(res_scal_respreco_prev);
  optim.d(:)=-1.*optim.residual_preco(:)+beta*optim.d(:);                  
  %Update iteration counter 
  optim.cpt_iter_CG=optim.cpt_iter_CG+1 ; 
  %-----------------------------------------------------%
  % if the debug option is activated, compute the       %
  % quadratic function minimized during the conjugate   %
  % gradient process (check is this value decresae      %
  % throughout the CG iterations )                      %
  %-----------------------------------------------------%  
  if(optim.debug)
     optim.qkm1_CG=optim.qk_CG;
     
     mgrad=-1.*grad(:)
     grad_term= scalL2(n,optim.descent,mgrad);
     descent_scal_Hd= scalL2(n,optim.descent_prev,optim.Hd);
     optim.hessian_term=optim.hessian_term+(optim.alpha_CG.^2)*optim.dHd+...
          2.*optim.alpha_CG*descent_scal_Hd;
     optim.qk_CG=-grad_term+0.5*optim.hessian_term;
  end
  %-----------------------------------------------------%
  % Check if the Eisenstat stopping critertion is       %
  % satisfied                                           %
  %-----------------------------------------------------% 
  optim.norm_residual= normL2(n,optim.residual);
  optim.conv_CG=((optim.norm_residual<=(optim.eta*optim.norm_grad))||(optim.cpt_iter_CG>=optim.niter_max_CG)) ;       
  %-----------------------------------------------------%
  % Print information on the current CG iteration       %
  %-----------------------------------------------------%
  % print_info_PTRN(optim,0e0,FLAG)     
  
end % descent_PTRN2
