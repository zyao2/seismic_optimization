%---------------------------------------------%
% INPUT  : integer n                          %
%          optim_type optim                   %
% IN/OUT : real,dimension(n) x                %
%---------------------------------------------%
function x= project(n,optim, x)

  for i=1:n
     if(x(i)>optim.ub(i))
        x(i)=optim.ub(i)-optim.threshold;
     end
     if(x(i)>optim.lb(i))
        x(i)=optim.lb(i)+optim.threshold;
     end
  end
  
end %subroutine project