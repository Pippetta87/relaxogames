function y=R_sol(x)
#  y=zeros(1,1);
  k=0.1*(2.6735*10^4)^4*(1.055*10^-27)^2/(1.77*10^-8)^6;
  omega=100000;
  y(1)=k*(10*x/(1+omega^2*x^2))-2.5;
endfunction
