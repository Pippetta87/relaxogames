%Create x and y values
xdata=omegaexp;
ydata=rexp;
w0=linspace(min(xdata),max(xdata),1000);
xdatal=log(omegaexp);
ydatal=log(rexp);
%Define the exponential model
%cgs: gauss=10^4T mu0=1 J=10^7erg
my_fun = @(p) sum((ydata-3/10*((gammaH*10^-4)^2*(10^7*hbarSI)/(p(1)^3*(10^-8)^3))^2*(p(2)./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2/ydata);
my_funl = @(p) sum((ydatal-log(3/10*(gammaH*10^-4)^2*10^7*hbarSI)^2+6*log(p(1)*(10^-8)^6)-log(p(2))+log(1./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2);
k
my_funk = @(p) sum((ydata-p(3)*(p(2)./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2/ydata);

%Function to evaluate the model. In this case we need to provide the error
%funtion: error sum of squares
fun_eval = @(p, w0) 3/10*((gammaH*10^-4)^2*hbarSI*10^7/(p(1)^3*(10^-8)^3))^2*(p(2)./(1+(w0*p(2)).^2)+2*p(2)./(1+4*(w0*p(2)).^2));
fun_evall = @(p, w0) log(3/10)+2*log((10^-4*gammaH)^2*hbarSI*10^7)-6*log(p(1)*(10^-8)^6)+log(p(2))-log(1./(1+(w0*p(2)).^2)+2*p(2)./(1+4*(w0*p(2)).^2));
fun_evalk = @(p, w0) p(3)*(p(2)./(1+(w0*p(2)).^2)+2*p(2)./(1+4*(w0*p(2)).^2));
%3/10*(mu0/(4*pi)*gH^2*hbar/(p(1)^3))^2
%Inital guesses
p0 = [1.5, 10^-7];
 
%Nonlinear fit
result = fminsearch(my_fun, p0)
%,'MaxFunEvals',10^4,'MaxIter',10^4)
3/10*((10^-4*gammaH)^2*hbarSI*10^7/(result(1)^3))^2
3/10*((10^-4*gammaH)^2*hbarSI*10^7/((10^-8)^3))^2
plot(xdata, ydata, 'k.')
hold on
plot(xdata, fun_eval(result, xdata), 'r')
xlabel('Time (s)')
ylabel('Signal')
