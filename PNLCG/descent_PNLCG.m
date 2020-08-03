% Computation of the descent direction given by the 
% preconditioned nonlinear conjugate gradient algorithm 
% of Dai and Yuan 
% Y. DAI AND Y. YUAN, A nonlinear conjugate gradient  %
% method with a strong global convergence property,   %
% SIAM Journal on Optimization, 10 (1999), pp. 177–182%
%                                                     %
% See also Nocedal, Numerical optimization,           %
% 2nd edition p.132                                   %
%-----------------------------------------------------%
% INPUT  : integer :: n (dimension)                   % 
%          real,dimension(n) :: grad,grad_preco       %
% INPUT/OUTPUT : optim_typ optim (data structure)     %
%-----------------------------------------------------%
function optim= descent_PNLCG(n,grad,grad_preco,optim)

  
  %------------------------------------------------------------%
  % Storing old descent direction                              %
  %------------------------------------------------------------%
  optim.descent_prev(:)=optim.descent(:);
  
  %------------------------------------------------------------%
  % Computation of beta                                        %
  %------------------------------------------------------------%  % 
  gkpgk= scalL2(n,grad,grad_preco);
  sk=zeros(n,1);
  sk(:)=grad(:)-optim.grad_prev(:);
  skpk= scalL2(n,sk,optim.descent_prev);
  beta=gkpgk/skpk;
  
  %------------------------------------------------------------%
  % Safeguard (may be useful in some cases)                    %
  %------------------------------------------------------------% 
  if((beta>=1e5)||(beta<=-1e5))
     beta=0.;
  end
  
  %------------------------------------------------------------%
  % Computation of the descent direction                       %
  %------------------------------------------------------------% 
  optim.descent(:)=-1.*grad_preco(:)+beta*optim.descent_prev(:);
  
  %------------------------------------------------------------%
  % Deallocation                                               % 
  %------------------------------------------------------------%
  clear sk;

end % descent_PNLCG
%