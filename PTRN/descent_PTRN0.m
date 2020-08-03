%-----------------------------------------------------%
% descent_PTRN0
%-----------------------------------------------------%
% INPUT : integer n dimension                         %
%         real,dimension(n) grad current gradient     %
%         real,dimension(n) grad_preco preconditioned %
%                           current gradient          %
%         character*4 FLAG communication flag         %
% INPUT/OUTPUT : optim_type optim (data structure)    %
%-----------------------------------------------------%
function optim=descent_PTRN0(n,grad,grad_preco,optim,FLAG)
false=0;
  %-----------------------------------------------------%
  % Initialization of the conjugate gradient process    %
  %-----------------------------------------------------%
  optim.residual(:)=grad(:);
  optim.residual_preco(:)=grad_preco(:);
  optim.d(:)=-1.*optim.residual_preco(:);
  optim.Hd(:)=0.;
  optim.descent(:)=0. ;    
  optim.qk_CG=0.;
  optim.hessian_term=0.;     
  optim.norm_residual=normL2(n,optim.residual);
  optim.conv_CG=false;    
  optim.cpt_iter_CG=0;
  %print_info_PTRN(optim,0e0,FLAG)     

 end % descent_PTRN0

