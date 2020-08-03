%*********************************************%
%*    SEISCOPE OPTIMIZATION TOOLBOX          *%
%*********************************************%
% This routine is used to deallocate all the  %
% arrays that have been used during the       %
% PSTD optimization                           %
%---------------------------------------------%
% INPUT/OUT:  optim_type optim                %
%---------------------------------------------%
function optim=finalize_PTRN(optim)
 clear optim.xk;
  clear optim.grad;
  clear optim.descent;
  clear optim.descent_prev;
  clear optim.residual;
  clear optim.residual_preco;
  clear optim.d;
  clear optim.Hd;
  clear optim.eisenvect; 
 
  
end % finalize_PTRN