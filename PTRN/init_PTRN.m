%------------------------------------------------------
% INPUT : integer n dimension                         %
%         real,dimension(n) x current point           %
%         real,dimension(n) grad current gradient     %
%         character*4 FLAG communication flag         %
% INPUT/OUTPUT : optim_type optim (data structure)    %
%-----------------------------------------------------%
function optim= init_PTRN(n,x,fcost,grad,optim)
true=1;
  %---------------------------------------%
  % set counters                          %
  %---------------------------------------%
  optim.cpt_iter=0;
  optim.f0=fcost;
  optim.nfwd_pb=0;
  optim.nhess=0;
  optim.eta=0.9 ; 
  
  %---------------------------------------%
  % initialize linesearch parameters      %
  % by default, the max number of         %
  % linesearch iteration is set to 20     %
  % and the initial steplength is set to 1%
  %---------------------------------------% 
  optim.m1=1e-4;  % Wolfe conditions parameter 1 (Nocedal value)
  optim.m2=0.9;   % Wolfe conditions parameter 2 (Nocedal value)
  optim.mult_factor=10; %Bracketting parameter (Gilbert value)
  optim.fk=fcost;
  optim.nls_max=20; % max number of linesearch
  optim.cpt_ls=0;
  optim.first_ls=true;
  optim.alpha=1.;   % first value for the linesearch steplength

  %---------------------------------------%
  % memory allocations                    %
  %---------------------------------------%
  optim.xk=zeros(n,1);
  optim.xk(:)=x(:);     
  optim.grad=zeros(n,1);
  optim.grad(:)=grad(:);
  optim.descent=zeros(n,1);
  optim.descent_prev=zeros(n,1);
  optim.residual=zeros(n,1);
  optim.residual_preco=zeros(n,1);
  optim.d=zeros(n,1);
  optim.Hd=zeros(n,1);
  optim.eisenvect=zeros(n,1);

  %---------------------------------------%
  % norm of the first gradient            %
  %---------------------------------------%
  optim.norm_grad=normL2(n,grad);
   
  
end % init_PTRN
