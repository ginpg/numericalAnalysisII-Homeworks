#! octave-interpreter-name -qf

# Giselt Parra, 26609640

tolerance     = 1e-6;     
epsilon       = 1e-14;     
maxIterations = 40;      


function  [xk1,i,l_error] = newton(f,d,xk,MaxIter,error,e)
  l_error = [];
  for i = 0:MaxIter;
%    if (feval(d,xk) <= e)
%      break
%    end
    pk = (-1) * feval(f,xk)/feval(d,xk);
    xk1 = xk + pk;
    l_error = cat(2,abs(feval(f,xk1)),l_error);
    if abs(feval(f,xk1)) <= error;
      break  
    end
    xk = xk1;
  end
end


function [xk1,i,l_error] = secant(f,xk,xk1,MaxIter,error,e)
  l_error = [];
  for i = 0:MaxIter
    xk = xk;
    xk1 = xk1;
%    if abs(xk - xk1) <= e
%      break  
%    end
    ak = (feval(f,xk) - feval(f,xk1))/ (xk - xk1); 
    if abs(ak) <= error;
      break  
    end      
    pk = (-1) *(feval(f,xk1)/ak);
    xk = xk1;
    xk1 = xk1 + pk;
    l_error = cat(2,abs(feval(f,xk1)),l_error);
    if abs(feval(f,xk1)) <= error;
      break  
    end
  end
end

function [new_point,i,l_error] = best_secant(f,xk,xk1,MaxIter,error,e)
  l_error = [];
  delta  = 1e-6; 
  alpha  = 1e-6; 
  a = xk;
  b = xk1;
  for i = 0:MaxIter
    % Escogencia del más cercano a la raiz
    fa = feval(f,a);
    fb = feval(f,b);
    if fa*fb < 0
      if abs(fa) < abs(fb)
        select = a;
        if a < b
          delta = abs(a-b)/100;
        else
          delta = (-1) * abs(a-b)/100;
        end
      else 
        select = b;
        if a > b
          delta = abs(a-b)/100;
         else
          delta = (-1) * abs(a-b)/100;
        end
      end
      
     elseif abs(fa) < abs(fb)
        select = a;
        if a < b
          delta = (-1) * abs(a-b)/100;
        else
          delta =  abs(a-b)/100;
        end
      else 
        select = b;
        if a > b
          delta = (-1) * abs(a-b)/100;
         else
          delta =  abs(a-b)/100;
        end
      end
    new_point = select;
  
    fsel = feval(f,select);
	  ak = (-1) * (select*delta*fsel)/(feval(f,select+delta*select) - fsel);
    
    if abs(ak) <= error
      break;
    end      
    
    new_point = select + ak;
    l_error = cat(2,abs(feval(f,new_point)),l_error);

    if abs(feval(f,new_point)) <= error
      break;
    end
	  a = new_point;
	  b = select;
  end
end


function  h = nro_iter(f,fprime,a,b,maxIterations,tolerance,epsilon,n);
  
  #Descomentar los disp(root_x) si desea verificar que el metodo halla la raiz de la funcion
  
  [root_sec,l7,e_sec] = secant(f,a,b,maxIterations,tolerance,epsilon);
    disp("Iteraciones necesitadas con Secante");
    disp(l7);  
    #disp(root_sec); 
    
   [root_new,l10,e_new] = newton(f,fprime,a,maxIterations,tolerance,epsilon);
    disp("Iteraciones necesitadas con newton");
    disp(l10); 
    #disp(root_new);
    
  [root_bestsec,l9,e_bestsec] = best_secant(f,a,b,maxIterations,tolerance,epsilon);
    disp("Iteraciones necesitadas con secante (BEST VERSION)");
    disp(l9); 
    #disp(root_bestsec);
    

  max_i = max(max(l7,l10),l9)+1;
  
  e_sec = flip(e_sec);
  e_new = flip(e_new);
  e_bestsec = flip(e_bestsec);

  e_sec(end+1: max_i) = 0;
  e_new(end+1: max_i) = 0;
  e_bestsec(end+1: max_i) = 0;


  a = [e_sec;e_new;e_bestsec];

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#DESCOMENTAR BLOQUE SI DESEA VER LA GRAFICA DEL RESIDUO DE LAS TRES EJECUCIONES EN UNA SOLA FIGURA
#DESCOMENTAR LINEA 230 ("f1 = figure...")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  subplot(3,3,n);
%  h = semilogy(a','-o','LineWidth',2);
%  set(gca, 'box', 'off');
%  xlabel('k');
%  ylabel('Error');
%  %xlim([1 15]);
%  
%  if (n == 1)
%    title({'x^3 + 2*x^2 - x + 4','a=x0=-50, b=x1=-20','root = -2.8455',' ','Gráfica del error con escala logarítmica','de base 10 en el eje y'});
%  elseif (n == 2)
%    title({'x^3 + 2*x^2 - x + 4','a=x0=-3, b=x1=-2','root = -2.8455',' ','Gráfica del error con escala logarítmica','de base 10 en el eje y'});
%  else
%    title({'x^2 - x - 1','a=x0=1,  b=x1=3','root = -2.8455',' ','Gráfica del error con escala logarítmica','de base 10 en el eje y'});
%  end
%  
%  subplot(3,3,n+3);
%  h1 = semilogx(a','-o','LineWidth',2);
%  set(gca, 'box', 'off');
%  xlabel('k');
%  ylabel('Error');
%  xlim([1 max_i]);
%  title({'Gráfica del error con escala logarítmica','de base 10 en el eje x'});
%  
%   
%  if (n == 3)
%    hL = subplot(3,3,8);
%    poshL = get(hL,'position');     
%    lgd = legend(hL,h,'secante','newton','secante mejorado');
%    set(lgd,'position',poshL);      % Adjusting legend's position
%    axis(hL,'off');                 % Turning its axis off
%    
%  end



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#DESCOMENTAR BLOQUE SI DESEA VER LA GRAFICA DEL RESIDUO CON UNA FIGURA POR PRUEBA

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  figure('Position', [600 600 900 300]);
%
%  subplot(1,3,2);
%  h1 = semilogy(a','-o','LineWidth',2);
%  xlabel('k');
%  ylabel('Error');
%  title({'Gráfica del error con escala logarítmica','de base 10 en el eje y'});
%  set(gca, 'box', 'off');
%  xlim([1 max_i]);
%  
%
%  subplot(1,3,3);
%  h2 = semilogx(a','-o','LineWidth',2);
%  xlabel('k');
%  ylabel('Error');
%  title({'Gráfica del error con escala logarítmica','de base 10 en el eje x'});
%  set(gca, 'box', 'off');
%  xlim([1 max_i]);
%  
%  hL = subplot(1,3,1);
%  poshL = get(hL,'position');    
%  lgd = legend(hL,[h1;h2],'secante','newton','secante mejorado');
%  set(lgd,'position',poshL);      % Adjusting legend's position
%  axis(hL,'off');                 % Turning its axis off
%
%  if (n == 1)
%    title({'x^3 + 2*x^2 - x + 4','a=x0=-50, b=x1=-20'});
%  elseif (n == 2)
%    title({'x^3 + 2*x^2 - x + 4','a=x0=-3, b=x1=-2'});
%  else
%    title({'x^2 - x - 1','a=x0=1,  b=x1=3'});
%  end
%  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%f1 = figure('Position', [100 500 1000 700],'Color', [.8 .8 .8]);



disp("\nDATA 0");
a             = -50;
b             = -20;
f             = @(x) x^3 + 2*x^2 - x + 4;
fprime        = @(x) 3*x^2 + 4*x - 1;
g             = @(x) nthroot((-2*x^2 + x - 1),3);

disp(f);
fprintf('a=x0=%d, b=x1=%d\n',a,b);
nro_iter(f,fprime,a,b,maxIterations,tolerance,epsilon,1);


disp("\nDATA 1");
a             = -3;
b             = -2;
f             = @(x) x^3 + 2*x^2 - x + 4;
fprime        = @(x) 3*x^2 + 4*x - 1;
g             = @(x) nthroot((-2*x^2 + x - 1),3);

disp(f);
fprintf('a=x0=%d, b=x1=%d\n',a,b);
nro_iter(f,fprime,a,b,maxIterations,tolerance,epsilon,2);


disp("\nDATA 2");
a             = 1;
b             = 3;
f             = @(x) x^2 - x - 1;
fprime        = @(x) 2*x - 1;
g             = @(x) 1 + 1/x ;

disp(f);
fprintf('a=x0=%d, b=x1=%d\n',a,b);
nro_iter(f,fprime,a,b,maxIterations,tolerance,epsilon,3);



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
