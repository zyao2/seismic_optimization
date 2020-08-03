%---------------------------------------------%
% INPUT/OUT:  optim_type optim                %
%---------------------------------------------%
function  optim=finalize_TRN(optim)

  clear optim.xk;
  clear optim.grad;
  clear optim.descent;
  clear optim.descent_prev;
  clear optim.residual;
  clear optim.d;
  clear optim.Hd;
  clear optim.eisenvect; 
  
end %subroutine finalize_TRN