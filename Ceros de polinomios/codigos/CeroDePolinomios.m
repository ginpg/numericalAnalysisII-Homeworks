#! octave-interpreter-name -qf

# Giselt Parra, 26.609.640

tolerance     = 1e-6;    

%Deflacion
function [y,new_coef] = horner_deflation(coef,r);
  new_coef = [];
  new_coef = cat(2,coef(1),new_coef);
  for i = 2:length(coef);
    new_coef = cat(2,coef(i) + new_coef(1)*r,new_coef);
  end
  y = new_coef(1);
  new_coef = flip(new_coef);
end


function root = solver(a,b,c,x2);
  p1 = sqrt(b^2 - 4*a*c);
  root1 = (-1)*b + p1;
  root2 = (-1)*b - p1;
  root1 = root1/(2*a);
  root2 = root2/(2*a);
  dif = abs(root1-x2);
  if dif < abs(root2-x2);
    root = root1;
  else
    root = root2;
  end
end


function root = lagrange_inte(x0,x1,x2,fx0,fx1,fx2);
  %Denominadores
  d0 = (x0-x1)*(x0-x2);
  d1 = (x1-x0)*(x1-x2);
  d2 = (x2-x0)*(x2-x1);
  
  m1 =  fx0/d0;
  m2 =  fx1/d1;
  m3 =  fx2/d2;
  
  %Coeficientes del polinomio interpolante
  a = m1 + m2 + m3;
  b = (-1) * (m1*(x1+x2) + m2*(x0+x2) + m3*(x0+x1));
  c = m1*x1*x2 + m2*x0*x2 + m3*x0*x1;

  root = solver(a,b,c,x2);
end
   


function roots = principal(poly,a,b,c,roots);
  while (length(poly) > 2)

    fa = horner_deflation(poly, a)(1);
    fb = horner_deflation(poly, b)(1);
    fc = horner_deflation(poly, c)(1);
    new_point = lagrange_inte(a,b,c,fa,fb,fc); 
    
    if abs(c-new_point) <= 1e-6;  
      [dif, poly] = horner_deflation(poly, new_point);
      poly(end) = [];
      roots = cat(2,new_point,roots);
      if !isreal(new_point);
        [dif, poly] = horner_deflation(poly, conj(new_point));
        poly(end) = [];
        roots = cat(2,conj(new_point),roots);
        %disp("Es raiz imaginaria");
        continue;
      end
      %disp("Es raiz real");
    else
      a = b;
      b = c;
      c = new_point;
    end
  end
  roots = cat(2,poly(2),roots);
end


% 2.6906, −0.3453 ± 1.31876i
disp("p(x) = x^3 - 2x^2 - 5");
p1 = [1,-2,0,-5]
[roots] = principal([1,-2,0,-5],-1,0,1,[])


% 5.8210, 1.6872, −0.5090
disp("p(x) = x^3 - 7x^2 + 6x + 5");
p2 = [1,-7,6,5]
[roots] = principal([1,-7,6,5],0,1,2,[])


x=-5:.1:5; % polynomial function
plot(x,polyval(p1,x))
grid on

x=-5:.1:5; % polynomial function
plot(x,polyval(p1,x))
grid on


figure('Position', [600 600 900 300]);

plot(x,polyval(p1,x))
xlabel('x');
ylabel('p_3(x)');
title({'p(x) = x^3 - 2x^2 - 5'});
set(gca, 'box', 'off');




