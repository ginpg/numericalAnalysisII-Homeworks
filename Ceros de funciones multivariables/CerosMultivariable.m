#! octave-interpreter-name -qf

# Giselt Parra, 26609640

tolerance     = 1e-6;      
maxIterations = 300;      


%F(x)
function m = F1(x,n)
  if n == 1
    m = [x(1)^2 - 2*x(2) + 1; 2*x(1) + x(2)^2 - 3];
  elseif n == 2
    m = [x(1)^2 + x(1)*x(2)^3 - 9  ;3*x(1)*x(2) - x(2)^3 - 4];
  elseif n == 3
    a = x(1)^2 + x(2)^2 -1;
    b = sin(pi*0.5*x(1))+ x(2)^3;
    m = [a; b];
  end

end

%J(x)
function m = A1(x,n)
  if n == 1 
    m = [2*x(1)  (-2); 2  2*x(2)];
  elseif n ==2
    m = [2*x(1)+x(2)^3  3*x(1)*x(2)^2; 3*x(2)   3*(x(1)-x(2)^2)];
  elseif n ==3
    m = [2*x(1)  2*x(2); pi*cos(pi*x(1)/2)/2     3*x(2)^2];    
  end
end

%Convergencia: q-cuadratica
function [i,x1,l_error] = Newton(x, MaxIter, error,n)
  l_error = [];
  newF = F1(x',n);
  for i = 0:MaxIter
    F = newF;
    A = A1(x',n);
    p = -A\F;
    x1 = x + p;
    newF = F1(x1,n);
    l_error = cat(2,norm(newF),l_error);
    l_error(1);
    if(l_error(1) < error);
      break
    end
    x = x1;
  end
end

%Convergencia: q-lineal
function [i,x1,l_error] = NewtonModificado(x, MaxIter, error,n)
  l_error = [];
  newF = F1(x',n);
  A = A1(x',n);
  for i = 0:MaxIter
    F = newF;
    p = -A\F;
    x1 = x + p;
    newF = F1(x1,n);
    l_error = cat(2,norm(newF),l_error);
    l_error(1);
    if(l_error(1) < error);
      break
    end
    x = x1;
  end
end

%function [i,x1,l_error] = DiferenciasFinitas(x, MaxIter, error,n)
%  l_error = [];
%  newF = F1(x',n);
%  hk = norm(newF);
%  A = (F1(abs(hk)*ones(1,length(x)),n) - newF)/hk
%  A = A*[1 2];
%  for i = 0:MaxIter
%    F = newF;
%    p = -A\F;
%    x1 = x + p;
%    newF = F1(x1,n);
%    l_error = cat(2,norm(newF),l_error);
%    if(l_error(1) < error);
%      break
%    end
%    x = x1;
%  end
%end



  function [i,x1,l_error,x_error] = Broyden(x, MaxIter, error,n)
  l_error = [];
  x_error = [];
  A = A1(x',n);
  newF = F1(x',n);
  for i = 0:MaxIter
    F = newF;
    p = -A\F;
    x1 = x + p;
    newF = F1(x1,n);
    x_error = cat(2,x1,x_error);
    l_error = cat(2,norm(newF),l_error);
    l_error(1);
    if(l_error(1) < error)
      break
    end
    y = newF - F;
    A = A +((y-A*p)*p')/(p'*p);
    x = x1;
  end
end



function [errorn,errorb,xerrorb] = zeros(x, maxIterations, tolerance,n)
  root_bro = [];
  root_new = [];
  [l1, rootsn, errorn]=  Newton(x, maxIterations, tolerance,n);
  [l11, rootsnm, errornm]=  NewtonModificado(x, maxIterations, tolerance,n);
  [l2, rootsb, errorb, xerrorb] = Broyden(x, maxIterations, tolerance,n);

  if n == 1
    F = "\nx² - 2y + 1 \n2x + y² - 3";
  elseif n == 2
    F = "\nx² + xy³ - 9 \n3xy - y³ - 4";
  elseif n == 3
    F = "\nx² + y² -1, sin(pi.x/2)+ y³";
  end
  
  F=F
  x= x'
  disp("Iteraciones necesitadas con Newton:");
  disp(l1);

  disp("Iteraciones necesitadas con Newton Modificado");
  disp(l11);
  
  disp("Iteraciones necesitadas con Broyden:");
  disp(l2);
  disp("\n\n");
  
  max_i = max(l1,l2)+1;
  errorb = flip(errorb);
  errorn = flip(errorn);
  errorb(end+1: max_i) = 0;
  xerrorb(end+1: max_i) = 0;
  errorn(end+1: max_i) = 0;
end



x = [2;2];
[e_new2,f_bro_error, x_bro_aprox]  = zeros(x, maxIterations, tolerance,1);

x = [1.5;1.5];
[e_new2,f_bro_error, x_bro_aprox]  = zeros(x, maxIterations, tolerance,1);


x = [1.5;1.5];
[e_new2,f_bro_error, x_bro_aprox] = zeros(x, maxIterations, tolerance,2);

x = [2;2];
[e_new2,f_bro_error, x_bro_aprox] = zeros(x, maxIterations, tolerance,2);

x = [1;1];
[e_new2,f_bro_error, x_bro_aprox] = zeros(x, maxIterations, tolerance,3);
