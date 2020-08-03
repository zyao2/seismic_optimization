%-----------------------------------------------------%
% INPUT : integer n, dimension of the problem         %
%       : real fcost                                  %
%       : real,dimension(n) x,grad,grad_preco         %
%                           initial guess, gradient,  %
%                           preconditioned gradient   %
% INPUT/OUTPUT : optim_type optim data structure      % 
%-----------------------------------------------------%
function optim= init_PSTD(n,x,fcost,grad,grad_preco,optim)

  %---------------------------------------%
  % set counters                          %
  %---------------------------------------%
  optim.cpt_iter=0;
  optim.f0=fcost;
  optim.nfwd_pb=0;  

  %---------------------------------------%
  % initialize linesearch parameters      %
  % by default, the max number of         %
  % linesearch iteration is set to 20     %
  % and the initial steplength is set to 1%
  %---------------------------------------% 
  optim.m1=1e-4; % Wolfe conditions parameter 1 (Nocedal value)
  optim.m2=0.9;  % Wolfe conditions parameter 2 (Nocedal value)
  optim.mult_factor=10; % Bracketting parameter (Gilbert value)
  optim.fk=fcost;
  optim.nls_max=20; % max number of linesearch
  optim.cpt_ls=0;
  optim.first_ls=1;%.true.
  optim.alpha=1.; % first value for the linesearch steplength
  
  %---------------------------------------%
  % memory allocations                    %
  %---------------------------------------%
  optim.xk=zeros(n,1);
  optim.xk(:)=x(:);
  optim.grad=zeros(n,1);
  optim.grad(:)=grad(:);
  optim.descent=zeros(n,1);

  %---------------------------------------%
  % first descent direction               %
  %---------------------------------------%
  %optim.descent(:)=-1.*grad(:)
  optim.descent(:)=-1.*grad_preco(:);

end % init_PSTD
