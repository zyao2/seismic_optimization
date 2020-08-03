%-----------------------------------------------------%
% INPUT : integer n, dimension of the problem         %
%       : real fcost                                  %
%       : real,dimension(n) x,grad first iterate and  %
%                           gradient                  %
% INPUT/OUTPUT : optim_type optim data structure      % 
%-----------------------------------------------------%
function optim=init_LBFGS(n,l,x,fcost,grad,optim)

  %---------------------------------------%
  % set counters                          %
  %---------------------------------------%
  optim.cpt_iter=0;
  optim.f0=fcost;
  optim.nfwd_pb=0 ; 
  optim.sk=zeros(n,l);
  optim.yk=zeros(n,l);
  optim.sk(:,:)=0.;
  optim.yk(:,:)=0.;
  optim.cpt_lbfgs=1;
  
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
  optim.first_ls=1;%true;
  optim.alpha=1.; % first value for the linesearch steplength 
  optim.bound=0;
  %---------------------------------------%
  % memory allocations                    %
  %---------------------------------------%
  optim.xk=zeros(n,1);
  optim.xk(:)=x(:) ;  
  optim.grad=zeros(n,1);
  optim.grad(:)=grad(:);
  optim.descent=zeros(n,1);

  %---------------------------------------%
  % first descent direction               %
  %---------------------------------------%
  optim.descent(:)=-1.*grad(:) ; 
  optim=save_LBFGS(n,optim.cpt_lbfgs,l,x,grad,optim);
  
end % init_LBFGS
