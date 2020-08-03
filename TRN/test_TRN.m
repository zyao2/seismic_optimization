 %----------------------------------------------------%
 % parameter initialization                           %
 %----------------------------------------------------%
 addpath ../common
  n=2;                     % dimension
  FLAG='INIT';             % first flag
  optim.niter_max=100;     % maximum iteration number
  optim.conv=1e-8 ;        % tolerance for the stopping criterion
  optim.print_flag=1 ;     % print info in output files 
  optim.debug=0;     % level of details for output files
  optim.niter_max_CG=5;    % maximum number of inner conjugate gradient 
                          % iterations
  %----------------------------------------------------%
  % intial guess                                       %
  %----------------------------------------------------%
  x=zeros(n,1);
  grad=zeros(n,1);
  x(1)=2.;
  x(2)=1.5;

  %----------------------------------------------------%
  % computation of the cost and gradient associated    %
  % with the initial guess                             %
  %----------------------------------------------------%
  [fcost,grad]=rosenbrock(x);
 optim.bound=0;
  
  %----------------------------------------------------%
  % optimization loop: while convergence not reached or%
  % linesearch not failed, iterate                     %
  %----------------------------------------------------%
  while (~strcmp(FLAG,'CONV') && ~strcmp(FLAG,'FAIL'))
     [x,optim,FLAG]=TRN(n,x,fcost,grad,optim,FLAG);
     if(FLAG=='GRAD') 
        %----------------------------------------------------%
        % if FLAG is GRAD, then compute cost and gradient in %
        % fcost, grad                                        %
        %----------------------------------------------------%
        [fcost,grad]= rosenbrock(x);        
     elseif(strcmp(FLAG,'HESS'))       
        %----------------------------------------------------%
        % if FLAG is HESS, then multiply optim%d by the      %
        % Hessian operator and store the result in optim%Hd  %
        %----------------------------------------------------%
        optim.Hd= rosenbrock_hess(x,optim.d,optim.Hd);
     end
  end
  
  %Helpful console writings
 % write(*,*) 'END OF TEST'
 % write(*,*) 'FINAL iterate is : ', x(:)
 % write(*,*) 'See the convergence history in iterate_TRN.dat and iterate_TRN_CG.dat'