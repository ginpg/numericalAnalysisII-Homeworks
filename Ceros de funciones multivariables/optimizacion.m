#! octave-interpreter-name -qf

# Giselt Parra, 26609640

tolerance     = 1e-4;     
maxIterations = 40;      


%FUNCIONES CON SUS DERIVADAS PARCIALES
function m = G1(x,n)
  if n == 2
    %exp(x-1) + exp(1-y) + (x-y)^2
    m = [2*x(1) - 2*x(2) + exp(x(1)-1); -2*x(1) + 2*x(2) - exp(1-x(2))];
  elseif n == 3
    m = [x(1)^2 - 2*x(1) + x(2)^2 - x(3) + 1;
        x(1)*x(2)^2 - x(1) - 3*x(2) + x(2)*x(3) + 2;
        x(1)*x(3)^2 - 3*x(3) + x(2)*x(3)^2 + x(1)*x(2)];
  end
  
end


function m = H1(x,n)
  if n == 2
    m = [exp(x(1)-1)+2 (-2); -2  exp(1-x(2))+2];
  elseif n == 3
    m = [2*x(1) - 2   2*x(2)               (-1);
         x(2)^2-1     2*x(1)*x(2)-3+x(3)   x(2);
         x(3)^2+x(2)  x(3)^2+x(1)          2*x(1)*x(3)-3+2*x(2)*x(3)];
  end
end

function [i,x1,l_error] =  Newton(x, MaxIter, error,n)
  l_error = [];
  newG = G1(x',n);
  for i = 0:MaxIter
    G = newG;
    H = H1(x',n);
    p = -H\G;
    x1 = x + p;
    newG = G1(x1,n);
    l_error = cat(2,norm(newG),l_error);
    if(l_error(1) < error);
      break
    end
    x = x1;
  end
end


function [i,x1,l_error] = BFGS(x, MaxIter, error,n)
  l_error = [];
  H = H1(x',n);
  newG = G1(x',n);
  for i = 0:MaxIter
    G = newG;
    p = -H\G;
    x1 = x + p;
    newG = G1(x1,n);
    l_error = cat(2,norm(newG),l_error);
    if(l_error(1) < error)
      break
    end
    y = newG - G;
    H = H + (y*y')/(y'*p) + (newG*G')/(p'*G);
    x = x1;
  end
end

function [errorn,ebfgs] = zeros(x, maxIterations, tolerance,n)
  root_bfgs = [];
  root_new = [];
  [l1, rootsn, errorn]=  Newton(x, maxIterations, tolerance,n);
  [l2, rootsb, ebfgs] = BFGS(x, maxIterations, tolerance,n);

  if n == 2
    G = "\n[2x - 2y + e^(x-1)\n- 2x + 2y - exp(1-y)]";
  elseif n == 3
    G = "\n[x² - 2x + y² - z + 1 \nxy² - x - 3y + yz + 2\nxz² - 3z + yz² + xy]";
  end
  
  G 
  x
  disp("\nIteraciones necesitadas con Newton:");
  disp(l1);
%  disp("Root:");
%  disp(rootsn);


  disp("Iteraciones necesitadas con BFGS:");
  disp(l2);
%  disp("Root:");
%  disp(rootsb);
  
  max_i = max(l1,l2)+1;
  ebfgs = flip(ebfgs);
  errorn = flip(errorn);
  ebfgs(end+1: max_i) = 0;
  errorn(end+1: max_i) = 0;
end


x = [0;0];
zeros(x, maxIterations, tolerance,2);


x = [0;0;0];
zeros(x, maxIterations, tolerance,3);









