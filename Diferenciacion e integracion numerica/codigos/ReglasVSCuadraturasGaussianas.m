
#Giselt Parra, 26.609.640

xi  = 0.5;

function ans = f(x)
  ans = 4/(1+x^2);
end


function metodos(a,b,n)
  n = n
  h = (a+b)/n;
  s_r = f((a + (a+h))/2);
  s_t = f(a)/2;
  s_s = f(a);
  
  for i = 1:n-1;
    s_r += f(((a+i*h) + (a+((i+1)*h)))/2);
    s_t += f(a+i*h);
    if  rem(i,2)==0;
      s_s += 2*f(a+i*h);
    else
      s_s += 4*f(a+i*h);
    end
  end 

  s_t += f(b)/2;
  s_s += f(b);  
  
  s_r = s_r/n;
  s_t =  h*s_t;
  s_s = (h/3)*s_s;
  
  disp("Resultado por regla de rectangulo")
  disp(s_r)
  disp("Resultado por regla de trapecio")
  disp(s_t)
  disp("Resultado por regla de Simpson")
  disp(s_s)

  
  disp("Errores (rec/trap/Simp)")
  disp(abs(pi-s_r))
  disp(abs(pi-s_t))
  disp(abs(pi-s_s))
  disp("\n")
end

metodos(0,1,8);
metodos(0,1,32);
metodos(0,1,128);