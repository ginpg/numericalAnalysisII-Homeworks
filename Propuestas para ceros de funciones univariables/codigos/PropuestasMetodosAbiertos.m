#! octave-interpreter-name -qf

# Giselt Parra, 26609640

tolerance     = 1e-6;      
epsilon       = 1e-14;    
maxIterations = 40;      

#Data 1
x0            = -0.3;
x1            = 2;
a             = -0.5;
b             = 1.7;
f             = @(x) x^2 - x - 1;
fprime        = @(x) 2*x - 1;
g             = @(x) 1 + 1/x ;


function  [c,i,l_error] = bisection(f,a,b,MaxIter,error)
  l_error = [];
  for i = 0:MaxIter
    c = a + (b-a)/2;
    if (feval(f,a)*feval(f,c))>0 ;       #f(a)f(c)>0
      a = c;
    elseif feval(f,b)*feval(f,c)>0;      #f(b)f(c)>0
      b = c;
    end  
    l_error = cat(2,abs(b-a),l_error);
    if (abs(b-a) <= error) or (abs(feval(f,c)) <= error);
      break
    end
  end
end 

function  [c,i,l_error] = better_bisection(f,a,b,MaxIter,error)
  l_error = [];
  for i = 0:MaxIter
    c = a + (b-a)/3;
    d = c + (b-a)/3;
    
    if (feval(f,a)*feval(f,c))<=0    
      a = a;
      b = c;
    elseif (feval(f,b)*feval(f,d))<=0
      a = d;
      b = b;
    else
      a = c;
      b = d;
    end

    l_error = cat(2,abs(b-a),l_error); 
    abs(b-a);
    abs(b-a) <= error;
    
    feval(f,c);
    abs(feval(f,c)) <= error;
    
    feval(f,d);
    abs(feval(f,d)) <= error;
    
    if (abs(b-a) <= error);
      c = a + (b-a)/2;
      break   
    elseif (abs(feval(f,a)) <= error);
      c = a;
      break   
    elseif (abs(feval(f,b)) <= error);
      c = b;
      break
    end
  end
end 

function  [xk2,i,l_error] = regula_Falsi(f,xk,xk1,MaxIter,error,e)
  l_error = [];
  for i = 0:MaxIter
%    if abs((feval(f,xk1) - feval(f,xk))) <= error
%      break  
%    end   
    ak = (xk1 - xk)/(feval(f,xk1) - feval(f,xk));    
    pk = -1 * feval(f,xk) * ak;     
    xk2 = xk + pk;
    l_error = cat(2,abs(feval(f,xk2)),l_error);  
    if (abs(xk1 - xk2) <= error) or (abs(xk - xk2) <= error);
      break  
    elseif abs(feval(f,xk2)) <= error;
      break
    end
    if (feval(f,xk1)*feval(f,xk2)<0);
      xk = xk1;
      xk1 = xk2;
    else
      xk1 = xk2;
    end
  end
end

function  [a2,i,l_error] = better_regula_Falsi(f,xk,xk1,MaxIter,error,e)
  l_error = [];
  a = xk;
  b = xk1;
  for i = 0:MaxIter
    
    c = a + (b-a)/2;
    if (feval(f,a)*feval(f,c))>0 ;       
      a = c;
    elseif feval(f,b)*feval(f,c)>0;     
      b = c;
    end  
    a = a;
    b = b;
    
    ak = (b - a)/(feval(f,b) - feval(f,a));    
    pk = -1 * feval(f,a) * ak;    
    
    a2 = a + pk;
    
    l_error = cat(2,abs(feval(f,a2)),l_error);

    if (abs(b - a2) <= error) or (abs(a - a2) <= error); 
      a2;
      break  
    elseif abs(feval(f,a2)) <= error;
      a2;
      break
    end
    
    
    if (feval(f,b)*feval(f,a2)<0);
      a = b;
      b = a2;
    else
      b = a2;
    end
  end
end

function  [a2,i,l_error] = better_bisection_regula_Falsi(f,xk,xk1,MaxIter,error,e)
  l_error = [];
  a = xk;
  b = xk1;
  
  for i = 0:MaxIter
    c = a + (b-a)/3;
    d = c + (b-a)/3;
    
    if (feval(f,a)*feval(f,c))<=0    
      a = a;
      b = c;
    elseif (feval(f,b)*feval(f,d))<=0
      a = d;
      b = b;
    else
      a = c;
      b = d;
    end
     
    a = a;
    b = b;
    
    ak = (b - a)/(feval(f,b) - feval(f,a));    
    pk = -1 * feval(f,a) * ak;    
    
    a2 = a + pk;
    
    l_error = cat(2,abs(feval(f,a2)),l_error); 
   
    if (abs(b-a) <= error);
      %a2 = a + (b-a)/2;
      break   
    elseif (abs(feval(f,a)) <= error);
      a2 = a;
      break   
    elseif (abs(feval(f,b)) <= error);
      a2 = b;
      break
    elseif (abs(feval(f,a2)) <= error);  
      a2;
      break
    end


    if (abs(b - a2) <= error) or (abs(a - a2) <= error);   
      a2;
      break  
    end
    
    if (feval(f,b)*feval(f,a2)<0);
      a = b;
      b = a2;
    else
      b = a2;
    end
  end
end




function  h = nro_iter(f,a,b,maxIterations,tolerance,epsilon,n);
  
   [root_bi,l1,e_bi] = bisection(f,a,b,maxIterations,tolerance);
    disp("Iteraciones en metodo de la Biseccion:");
    disp(l1);

  [root_bbi,l3,e_bbi] = better_bisection(f,a,b,maxIterations,tolerance);
    disp("Iteraciones en metodo de la Biseccion mejorado:");
    disp(l3);
    
  [root_reg,l4,e_reg] = regula_Falsi(f,a,b,maxIterations,tolerance,epsilon);
    disp("Iteraciones en Regula Falsi:");
    disp(l4);
    
  [root_breg,l5,e_breg] = better_regula_Falsi(f,a,b,maxIterations,tolerance,epsilon);
    disp("Iteraciones en Regula Falsi mejorado:");
    disp(l5);
    
  [root_bireg,l6,e_bireg] = better_bisection_regula_Falsi(f,a,b,maxIterations,tolerance,epsilon);
    disp("Iteraciones necesitadas con Biseccion-Regula Falsi mejorado:");
    disp(l6);
    
    
  max_i = max(max(max(max(l1,l6),l3),l4),l5)+1;
  
  e_bi = flip(e_bi);
  e_bbi = flip(e_bbi);
  e_reg = flip(e_reg);
  e_breg = flip(e_breg);
  e_bireg = flip(e_bireg);
  
  e_bi(end+1: max_i) = 0;
  e_bbi(end+1: max_i) = 0;
  e_reg(end+1: max_i) = 0;
  e_breg(end+1: max_i) = 0;
  e_bireg(end+1: max_i) = 0;

  a = [e_bi;e_bbi;e_reg;e_breg;e_bireg];

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#DESCOMENTAR BLOQUE SI DESEA VER LA GRAFICA DEL ERROR DE LAS TRES EJECUCIONES EN UNA SOLA FIGURA
#DESCOMENTAR LINEA 286 ("f1 = figure...")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  subplot(3,3,n);
%  h = semilogy(a','-o','LineWidth',2);
%  set(gca, 'box', 'off');
%  xlabel('k');
%  ylabel('Error');
%  xlim([1 15]);
%  
%  if (n == 1)
%    title({'x^3 + 2*x^2 - x + 4','a=x0=0, b=x1=-5','root = -2.8455',' ','Gráfica del error con escala logarítmica','de base 10 en el eje y'});
%  elseif (n == 2)
%    title({'x^3 + 2*x^2 - x + 4','a=x0=-3, b=x1=-2','root = -2.8455',' ','Gráfica del error con escala logarítmica','de base 10 en el eje y'});
%  else
%    title({'x^2 - x - 1','a=x0=1,  b=x1=2','root = -2.8455',' ','Gráfica del error con escala logarítmica','de base 10 en el eje y'});
%  end
%  
%  subplot(3,3,n+3);
%  h = semilogx(a','-o','LineWidth',2);
%  set(gca, 'box', 'off');
%  xlabel('k');
%  ylabel('Error');
%  xlim([1 15]);
%  title({'Gráfica del error con escala logarítmica','de base 10 en el eje x'});
%  
%   
%  if (n == 3)
%    hL = subplot(3,3,8);
%    poshL = get(hL,'position');     
%    lgd = legend(hL,h,'biseccion','biseccion mejorada','regula falsi','regula falsi mejorada','biseccion-regula falsi mejorada');
%    set(lgd,'position',poshL);      % Adjusting legend's position
%    axis(hL,'off');                 % Turning its axis off
%    
%  end


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%f1 = figure('Position', [100 500 1000 700],'Color', [.8 .8 .8]);



disp("\nDATA 0");
a             = 0;
b             = -10;
f             = @(x) x^3 + 2*x^2 - x + 4;
fprime        = @(x) 3*x^2 + 4*x - 1;
g             = @(x) nthroot((-2*x^2 + x - 1),3);

disp(f);
fprintf('a=%d, b =%d\n',a,b);
nro_iter(f,a,b,maxIterations,tolerance,epsilon,1);


disp("\nDATA 1");
a             = -3;
b             = -2;
f             = @(x) x^3 + 2*x^2 - x + 4;
fprime        = @(x) 3*x^2 + 4*x - 1;
g             = @(x) nthroot((-2*x^2 + x - 1),3);

disp(f);
fprintf('a=%d, b =%d\n',a,b);
nro_iter(f,a,b,maxIterations,tolerance,epsilon,2);


disp("\nDATA 2");
a             = 1;
b             = 2;
f             = @(x) x^2 - x - 1;
fprime        = @(x) 2*x - 1;
g             = @(x) 1 + 1/x ;

disp(f);
fprintf('a=%d, b =%d\n',a,b);
nro_iter(f,a,b,maxIterations,tolerance,epsilon,3);



f0 = figure('Position', [100 500 9 3]);
a = 0;
a = 0;
a = 0;
close(f0);

f3 = figure('Position', [100 500 9 3]);
a = 0;
a = 0;
a = 0;
close(f3);
