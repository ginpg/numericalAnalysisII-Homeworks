
#Giselt Parra, 26.609.640

xi  = 0.5;


f   =  @(x) (4*x^10)/(4608*6) + (x^10)/(9*6) + (30*x^2)/6 + (30*x)/6 - 11;
d   =  @(x) (40*(x^9))/(4608*6)  + (10*(x^9))/(9*6) + (60*x)/6 + 30/6;


function xi = newton(f,d,xi)
  for i = 0:40;
    pk = (-1) * feval(f,xi)/feval(d,xi);
    xi = xi + pk;
    if abs(feval(f,xi)) <= 1e-7;
      break  
    end
  end
end

disp("La raiz obtenida es: ")
disp(newton(f,d,xi));