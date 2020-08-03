  %----------------------------------------------------%
  % parameter initialization                           %
  %----------------------------------------------------%
addpath ../common
  n=2;                     % dimension
  FLAG='INIT';             % first flag
  optim.niter_max=10000;   % maximum iteration number 
  optim.conv=1e-8;         % tolerance for the stopping criterion
  optim.print_flag=1;      % print info in output files 
  optim.debug=0;%.false.     % level of details for output files  
  optim.l=20;              % maximum number of stored pairs used for
                          % the l-BFGS approximation
  
  %----------------------------------------------------%
  % intial guess                                       %
  %----------------------------------------------------%
  
x=zeros(n,1);
grad=zeros(n,1);
grad_preco=zeros(n,1);
  x(1)=1.5;
  x(2)=1.5;
  
  %----------------------------------------------------%
  % computation of the cost and gradient associated    %
  % with the initial guess                             %
  %----------------------------------------------------%
  [fcost,grad]= rosenbrock(x);
  grad_preco(:)=1.*grad(:);
  
  %----------------------------------------------------%
  % optimization loop: while convergence not reached or%
  % linesearch not failed, iterate                     %
  %----------------------------------------------------%
  while (~strcmp(FLAG,'CONV') && ~strcmp(FLAG,'FAIL'))
     [x,optim,FLAG]= PLBFGS(n,x,fcost,grad,grad_preco,optim,FLAG);
     if(strcmp(FLAG,'GRAD'))        
        [fcost,grad]= rosenbrock(x);       
     end
     if(strcmp(FLAG,'PREC'))
        %apply preconditioning to optim.q_plb
        %if nothing is done, PLBFGS is equivalent to LBFGS
     end
  end
 