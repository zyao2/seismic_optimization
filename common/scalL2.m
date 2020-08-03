function scal_xy=scalL2(n,x,y)

  scal_xy=0.;
  for i=1:n
     scal_xy=scal_xy+x(i)*y(i);
  end
  
end %subroutine scalL2
%