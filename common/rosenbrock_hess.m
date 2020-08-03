%---------------------------------------------%
% The Rosenbrock function is defined by       %
% f(x1,x2)  = (1-x1)**2+100.*(x2-x1**2)**2    %
%---------------------------------------------%
% The gradient of the Rosenbrock function is  %
% dx1f(x1,x2)= -2(x1-1)-400x1*(x2-x1**2)      %
% dx2f(x1,x2)= 200(x2-x1**2)                  %
%---------------------------------------------%
% The Hessian operator of the Rosenbrok       %
% function is                                 %
% dx1x1f(x1,x2)=-2-400(x2-x1**2)+800x1**2     %
% dx1x2f(x1,x2)=-400x1                        %
% dx2x2f(x1,x2)=200                           %
%---------------------------------------------%


%---------------------------------------------%
%  The routine Rosenbrock_Hess returns        %
%  Hessian-vector product H(x)d in output Hd  %
%  for input parameters x and d               %
%  H is the Hessian matrix                    %
%  x=(x1,x2), d=(d1,d2) are two vector of R^2 %
%---------------------------------------------%
function [Hd]=rosenbrock_hess(x,d,Hd)

  Hd(1)=(1200.*x(1).^2-400.*x(2)+2.)*d(1)-400*x(1)*d(2);
  Hd(2)=-400.*x(1)*d(1)+200.*d(2);
  
end % rosenbrock_hess
