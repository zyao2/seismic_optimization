  %----------------------------------------------------%
  % parameter initialization                           %
  %----------------------------------------------------%
  addpath ../common
  n=2;                     % dimension
  FLAG='INIT';             % first flag
  optim.niter_max=10000;   % maximum iteration number 
  optim.conv=1e-8;         % tolerance for the stopping criterion
  optim.print_flag=1;      % print info in output files 
  optim.debug=0;%.false.;     % level of details for output files
  optim.l=20 ;             % maximum number of stored pairs used for
                          % the l-BFGS approximation
 optim.bound=0;
  
  %----------------------------------------------------%
  % intial guess                                       %
  %----------------------------------------------------%
 x=zeros(n,1);
 grad=zeros(n,1);
  x(1)=1.5;
  x(2)=1.5;

  %----------------------------------------------------%
  % computation of the cost and gradient associated    %
  % with the initial guess                             %
  %----------------------------------------------------%
  [fcost,grad]=rosenbrock(x);
  
  %----------------------------------------------------%
  % optimization loop: while convergence not reached or%
  % linesearch not failed, iterate                     %
  %----------------------------------------------------%
   while (~strcmp(FLAG,'CONV') && ~strcmp(FLAG,'FAIL'))
     [x,optim,FLAG]= LBFGS(n,x,fcost,grad,optim,FLAG);
     if(FLAG=='GRAD')         
        %compute cost and gradient at point x
        [fcost,grad]=rosenbrock(x);        
     end
  end
  
