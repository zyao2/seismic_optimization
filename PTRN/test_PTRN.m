addpath ../common
false=0;
  %----------------------------------------------------%
  % parameter initialization                           %
  %----------------------------------------------------%
  n=2;                     % dimension
  FLAG='INIT';             % first flag
  optim.niter_max=100;     % maximum iteration number
  optim.conv=1e-8 ;        % tolerance for the stopping criterion
  optim.print_flag=1;      % print info in output files 
  optim.debug=false;     % level of details for output files
  optim.niter_max_CG=5 ;   % maximum number of inner conjugate gradient 
                          % iterations
  
  %----------------------------------------------------%
  % intial guess                                       %
  %----------------------------------------------------%
  x=zeros(n,1);grad=zeros(n,1);grad_preco=zeros(n,1);
  x(1)=1.5;
  x(2)=1.5;
  
  %----------------------------------------------------%
  % computation of the cost and gradient associated    %
  % with the initial guess                             %
  %----------------------------------------------------%
  [fcost,grad]= rosenbrock(x);

  %----------------------------------------------------%
  % multiplication of the gradient by the              %
  % preconditioner (here preconditioner is identity:   %
  % no preconditioning)                                %
  %----------------------------------------------------%
  grad_preco(:)=grad(:) ; 
  
  %----------------------------------------------------%
  % optimization loop: while convergence not reached or%
  % linesearch not failed, iterate                     %
  %----------------------------------------------------%
  while (~strcmp(FLAG,'CONV') && ~strcmp(FLAG,'FAIL'))
     [x,optim,FLAG]= PTRN(n,x,fcost,grad,grad_preco,optim,FLAG);
     if(FLAG=='GRAD')     
        %----------------------------------------------------%
        % if FLAG is GRAD, then compute cost, gradient and   %
        % preconditioned gradient in fcost, grad, grad_preco %
        %----------------------------------------------------%
        [fcost,grad]= rosenbrock(x);        
        grad_preco(:)=grad(:);
     elseif(FLAG=='HESS')
        %----------------------------------------------------%
        % if FLAG is HESS, then multiply optim%d by the      %
        % Hessian operator and store the result in optim%Hd  %
        %----------------------------------------------------%
        optim.Hd= rosenbrock_hess(x,optim.d,optim.Hd);
     elseif(FLAG=='PREC') 
        %----------------------------------------------------%
        % if FLAG is PREC, then multiply optim%residual by   %
        % preconditioner and store the result in             %
        % optim%residual_preco                               %
        %----------------------------------------------------%
        optim.residual_preco(:)=optim.residual(:);
     end
  end
  

