% This linesearch enforces the Wolfe          %
% conditions: sufficient decrease             %
%             sufficient curvature            %
% The Wolfe conditions can be found in        %
% Nocedal, Numerical Optimization, 2nd        %
% edition, p.33                               %
%                                             %
% The linesearch method implemented here is   %
% based first on a bracketing strategy, then  %
% on a dichotomy algorithm. A full description%
% of this strategy can be found in            %
% Numerical Optimizationn Theoretical and     %
% Practical Aspects, J.F.Bonnans, J.C.Gilbert,%
% C. Lemaréchal, C.A. Sagastizábal,           %
% Springer-Verlag, Universitext               %
%                                             %
%---------------------------------------------%
% INPUT  : integer :: n (dimension)           %
%          real fcost (current cost)          %
%          real,dimension(n) grad             %
% OUTPUT : real,dimension(n) x                %
%          optim_typ optim (data structure)   %
%---------------------------------------------%

function [x,optim]= std_linesearch(n,x,fcost,grad,optim)
false=0;true=1;
  if(optim.first_ls)
     %---------------------------------------%
     % FIRST LINESEARCH: initialization step %
     %---------------------------------------%
     optim.fk=fcost;
     q0= scalL2(n,grad,optim.descent);
     optim.q0=q0;
     %set the search interval bounds to 0
     optim.alpha_L=0.;
     optim.alpha_R=0.;          
     optim.task='NEW_GRAD';
     optim.first_ls=false;
     optim.xk(:)=x(:);
     x(:)=optim.xk(:)+optim.alpha*optim.descent(:);
     %IF BOUNDS ACTIVATED, PROJECT x(:) TO THE FEASIBLE ENSEMBLE
     if(optim.bound==1)
        x=project(n,optim,x);
     end
     optim.cpt_ls=0;
  elseif( (optim.cpt_ls>=optim.nls_max) && (fcost<optim.fk))
     %-----------------------------------------------------------------------%
     % if the number of linesearch iteration outreaches the maximum allowed  %
     % but a decrease of the misfit is produced then accept the steplength   %
     %-----------------------------------------------------------------------%
     optim.task='NEW_STEP';        
     optim.first_ls=true;
     %Compute new x in the descent direction     
     x(:)=optim.xk(:)+optim.alpha*optim.descent(:);     
     %IF BOUNDS ACTIVATED, PROJECT x(:) INTO TO THE FEASIBLE ENSEMBLE
     if(optim.bound==1)
        x= project(n,optim,x);
     end
  elseif(optim.cpt_ls>=optim.nls_max) 
     %-----------------------------------------------------------------------%
     % if the number of linesearch iteration outreaches the maximum allowed  %
     % without decreasing the misfit then the linesearch has failed          %
     %-----------------------------------------------------------------------%
     optim.task='FAILURE%';
  else
     %-----------------------------------------------------------------------%
     % If not initialization step and number of linesearch iteration ok      %
     % then perform one linesearch iteration                                 %
     %-----------------------------------------------------------------------%
     optim.q=scalL2(n,grad,optim.descent);
     if( (fcost<=(optim.fk+(optim.m1*optim.alpha*optim.q0)))&&...
          (optim.q>=(optim.m2*optim.q0)) )
        %--------------------------------------------------------------------%
        % First test if the Wolfe conditions are satisfied with              %     
        % current steplength, if this is the case, linesearch                % 
        % ends here with success                                             %
        %--------------------------------------------------------------------%
        optim.task='NEW_STEP';
        optim.first_ls=true;

     elseif (fcost>(optim.fk+(optim.m1*optim.alpha*optim.q0)))
        %--------------------------------------------------------------------%
        % If the first condition is not satisfied then shrink the            %
        % search interval                                                    %
        %--------------------------------------------------------------------%
   
        optim.alpha_R=optim.alpha;
        new_alpha=(optim.alpha_L+optim.alpha_R)/2.;
        optim.alpha=new_alpha;        
        optim.task='NEW_GRAD';
        optim.cpt_ls=optim.cpt_ls+1;
     elseif( (fcost<= (optim.fk+(optim.m1*optim.alpha*optim.q0))) && ...
          (optim.q<(optim.m2*optim.q0) ) ) 
        %--------------------------------------------------------------------%
        % If the second condition is not satisfied then shrink the           %
        % search interval unless the right bound of the search interval      %
        % as not yet been defined                                            %
        %--------------------------------------------------------------------%
   
        optim.alpha_L=optim.alpha;
        if(optim.alpha_R~=0.) 
           new_alpha=(optim.alpha_L+optim.alpha_R)/2.;
        else
           new_alpha=optim.mult_factor*optim.alpha;
        end
        optim.alpha=new_alpha;        
        optim.task='NEW_GRAD';
        optim.cpt_ls=optim.cpt_ls+1;                        
     end
     %Compute new x in the descent direction
     x(:)=optim.xk(:)+optim.alpha*optim.descent(:);     
     %IF BOUNDS ACTIVATED, PROJECT x(:) TO THE FEASIBLE ENSEMBLE
     if(optim.bound==1)
        x= project(n,optim,x);
     end
  end
end %subroutine std_linesearch


