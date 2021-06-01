#! octave-interpreter-name -qf

# Giselt Parra, 26609640

tolerance     = 1e-7;      # 7 digit accuracy is desired / error
epsilon       = 1e-14;     # Do not divide by a number smaller than this  / use to compare denominator
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

function  [xk1,i,l_error] = fixed_Point(f,g,xk,MaxIter,error,e)
  l_error = [];
  for i = 0:MaxIter
    xk = xk;
    xk1 = feval(g,xk);
    l_error = cat(2,abs(xk - xk1),l_error);
    if (abs(feval(g,xk1) - xk1)) <= e  or (abs(xk - xk1) <= e);  #abs(xk - xk1) < e
      break  
    end
    xk = xk1;
  end
end

function [e_bi,e_new,e_sec,e_reg,e_fix,max_i] = cal(f,g,fprime,a,b,x0,x1,maxIterations,tolerance,epsilon)
  [root_bi,l1,e_bi] = bisection(f,a,b,maxIterations,tolerance);
  [root_new,l2,e_new] = newton(f,fprime,x1,maxIterations,tolerance,epsilon);
  [root_sec,l3,e_sec] = secant(f,x0,x1,maxIterations,tolerance,epsilon);
  [root_reg,l4,e_reg] = regula_Falsi(f,x0,x1,maxIterations,tolerance,epsilon);
  [root_fix,l5,e_fix] = fixed_Point(f,g,x1,maxIterations,tolerance,epsilon);
  
  max_i = max(max(max(max(l1,l2),l3),l4),l5)+1;
  e_bi = flip(e_bi);
  e_new = flip(e_new);
  e_sec = flip(e_sec);
  e_reg = flip(e_reg);
  e_fix = flip(e_fix);
  e_bi(end+1: max_i) = 0;
  e_new(end+1: max_i) = 0;
  e_sec(end+1: max_i) = 0;
  e_reg(end+1: max_i) = 0;
  e_fix(end+1: max_i) = 0;
  
end


[e1_bi,e1_new,e1_sec,e1_reg,e1_fix,max1_i] = cal(f,g,fprime,0,2,0,2,maxIterations,tolerance,epsilon);

#Data 2 
x0            = 1.5;
x1            = 1.4;
a             = 1;
b             = 2;
f             = @(x) x^2 - x - 1;
fprime        = @(x) 2*x - 1;
g             = @(x) 1 + 1/x ;
[e2_bi,e2_new,e2_sec,e2_reg,e2_fix,max2_i] = cal(f,g,fprime,a,b,x0,x1,maxIterations,tolerance,epsilon);

#Data 3
x0            = -2.1;
x1            =  -2.2;
a             = -3;
b             = -2;
f             = @(x) x^3 + 2*x^2 - x + 4;
fprime        = @(x) 3*x^2 + 4*x - 1;
g             = @(x) nthroot((-2*x^2 + x - 1),3);

[e3_bi,e3_new,e3_sec,e3_reg,e3_fix,max3_i] = cal(f,g,fprime,a,b,x0,x1,maxIterations,tolerance,epsilon);

#Data 4
x0            = 0.01;
x1            = 0.4;
a             = 0;
b             = 1;
f             = @(x) x^4 + 2*x^2 - x -1;
fprime        = @(x) 4*x^3 + 4*x - 1 ;
g             = @(x) (-2*x^2 + x + 1 )^(1/4);
[e4_bi,e4_new,e4_sec,e4_reg,e4_fix,max4_i] = cal(f,g,fprime,a,b,x0,x1,maxIterations,tolerance,epsilon);



a1 = [e1_bi;e1_new;e1_sec;e1_reg;e1_fix];
a2 = [e2_bi;e2_new;e2_sec;e2_reg;e2_fix];
a3 = [e3_bi;e3_new;e3_sec;e3_reg;e3_fix];
a4 = [e4_bi;e4_new;e4_sec;e4_reg;e4_fix];

# F1
f1 = figure('Renderer', 'painters', 'Position', [100 500 1100 700]);

x = -100:100;

subplot(2,3,1);
fplot(@(x) x.^2 - x - 1,[-100 100]);
set(gca, 'box', 'off');
title({'f_1(x) = x^2 - x - 1','x_* = 1.6180'});
xlabel('x');
ylabel('y');
xlim([-20 20]);
ylim([-20 100]);
%hLeg = legend('example');
%set(hLeg,'visible','off');
grid on;

subplot(2,3,2);
h1 = semilogx(a1','-o','LineWidth',2);
set(gca, 'box', 'off')
title({'Velocidad de convergencia','a=-0.5, b=1.7, x_0=-0.3'});
xlabel('k');
ylabel('Error');
#lgd = legend('biseccion','newton','secante','regula falsi','fixed point');
xlim([1 max1_i]);
%ylim([0 1]);

subplot(2,3,3);
semilogx(a1','-o','LineWidth',2);
title('Acercamiento velocidad de convergencia en f_1');
xlabel('k');
ylabel('Error');
xlim([1 5]);
%ylim([0 1]);

# F2


subplot(2,3,5);
h2 = semilogx(a2','-o','LineWidth',2);
set(gca, 'box', 'off')
title({'Velocidad de convergencia','a=1, b=2, x_0=1.5'});
xlabel('k');
ylabel('Error');
%lgd = legend('biseccion','newton','secante','regula falsi','fixed point');
xlim([1 max2_i]);
%ylim([0 1]);


subplot(2,3,6);
semilogx(a2','-o','LineWidth',2);
title('Acercamiento velocidad de convergencia en f_1');
xlabel('k');
ylabel('Error');
xlim([1 5]);
%ylim([0 1]);

hL = subplot(2,3,4);
poshL = get(hL,'position');     % Getting its position

lgd = legend(hL,[h1;h2],'biseccion','newton','secante','regula falsi','fixed point');
set(lgd,'position',poshL);      % Adjusting legend's position
axis(hL,'off');                 % Turning its axis off


# F3
f2 = figure('Renderer', 'painters', 'Position', [100 500 1100 700]);
subplot(2,3,1);
fplot(@(x) x.^3 + 2*x.^2 - x + 4,[-100 100]);
set(gca, 'box', 'off')
title({'f_3(x) = x^3 + 2*x^2 - x + 4','root = -2.8455'});
xlabel('x');
ylabel('y');
xlim([-10 10]);
ylim([-50 50]);
%hLeg = legend('example');
%set(hLeg,'visible','off');
grid on;

subplot(2,3,2);
semilogx(a3','-o','LineWidth',2);
set(gca, 'box', 'off')
title({'Velocidad de convergencia','a=-3, b=-2, x0=-2.1'});
xlabel('k');
ylabel('Error');
#lgd = legend('biseccion','newton','secante','regula falsi','fixed point');
xlim([2 40]);
ylim([0 1]);

subplot(2,3,3);
semilogx(a3','-o','LineWidth',2);
title('Acercamiento velocidad de convergencia en f_2');
xlabel('k');
ylabel('Error');
xlim([1 11]);
ylim([0 0.45]);



# F4
subplot(2,3,4);
fplot(@(x) x.^4 + 2*x.^2 - x -1,[-100 100]);
set(gca, 'box', 'off')
title({'f_3(x) = x^4 + 2*x^2 - x -1','root = 0.82511'});
xlabel('x');
ylabel('y');
xlim([-10 10]);
ylim([-20 100]);
%hLeg = legend('example');
%set(hLeg,'visible','off');
grid on;

subplot(2,3,5);
semilogx(a4','-o','LineWidth',2);
set(gca, 'box', 'off')
title({'Velocidad de convergencia','a=0, b=2, x_0=0.6'});
xlabel('k');
ylabel('Error');
lgd = legend('biseccion','newton','secante','regula falsi','fixed point');
xlim([1 max4_i]);
ylim([0 10]);

subplot(2,3,6);
semilogx(a4','-o','LineWidth',2);
title('Acercamiento velocidad de convergencia en f_3');
xlabel('k');
ylabel('Error');
xlim([1 40]);
ylim([0 3]);

