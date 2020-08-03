%---------------------------------------------%
% INPUT:  optim_type optim                    %
% OUTPUT: logical test_conv                   %
%---------------------------------------------%
function test_conv=std_test_conv(optim,fcost)

 test_conv=((fcost/optim.f0<optim.conv) ||(optim.cpt_iter>optim.niter_max));
  
end %subroutine std_test_conv