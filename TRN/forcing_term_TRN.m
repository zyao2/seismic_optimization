%------------------------------------------------------
% INPUT : integer n dimension                         %
%         real,dimension(n) grad current gradient     %
%         character*4 FLAG communication flag         %
% INPUT/OUTPUT : optim_type optim (data structure)    %
%-----------------------------------------------------%
function optim=forcing_term_TRN(n,grad,optim)

  %-----------------------------------------------------%
  % Computation of the forcing term optim%eta following %
  % the formula                                         %
  %-----------------------------------------------------%
  eta_save=optim.eta ; 
  optim.eisenvect(:)=grad(:)-optim.residual(:);
  norm_eisenvect=normL2(n,optim.eisenvect);
  optim.eta=norm_eisenvect/optim.norm_grad_m1;
  
  %-----------------------------------------------------%
  % Additional safeguard if optim%eta is too large      %       
  %-----------------------------------------------------%
  eta_save_power=eta_save.^((1.+sqrt(5.))/2.);
  if(eta_save_power>0.1) 
     optim.eta=max(optim.eta,eta_save_power);
  end
  if(optim.eta>1.) 
     optim.eta=0.9;     
  end
  
end %subroutine forcing_term_TRN
