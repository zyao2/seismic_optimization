%---------------------------------------------%
% INPUT  : integer n                          %
%          real,dimension(n) x                %
% OUTPUT : real norm_x                        %
%---------------------------------------------%
function norm_x=normL2(n,x)
  norm_x=0.;
  for i=1:n
     norm_x=norm_x+x(i).^2;
  end
  norm_x=sqrt(norm_x);
  
%end %subroutine normL2
