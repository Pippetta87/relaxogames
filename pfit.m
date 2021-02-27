function [rHH,tc]=pfit(xdata,ydata,rHH0,tc0)
%Create x and y values
global gammaH;
global hbarSI;
%xdata=omegaexp;
%ydata=rexp;
%w0=linspace(min(xdata),max(xdata),1000);
xdatal=log(xdata);
ydatal=log(ydata);
%Define the exponential model
%cgs: gauss=10^4T mu0=1 J=10^7erg
lor2par = @(p) sum((ydata-3/10*((gammaH)^2*(hbarSI)/(p(1)^3*(10^-8)^3))^2*(p(2)./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2);
lorlog = @(p) sum((ydatal-log(3/10*(gammaH)^2*hbarSI)^2+6*log(p(1)*(10^-8)^6)-log(p(2))+log(1./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2);
lor1par = @(p) sum((ydata-p(3)*(p(2)./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2/ydata);
lor2parwou = @(p) sum((ydata-p(1)*(p(2)./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2);
%Function to evaluate the model. In this case we need to provide the error
%funtion: error sum of squares
lor2par_eval = @(p, w0) 3/10*((gammaH)^2*hbarSI/(p(1)^3*(10^-8)^3))^2*(p(2)./(1+(w0*p(2)).^2)+2*p(2)./(1+4*(w0*p(2)).^2));
lorlog_eval = @(p, w0) log(3/10)+2*log((gammaH)^2*hbarSI)-6*log(p(1)*(10^-8)^6)+log(p(2))-log(1./(1+(w0*p(2)).^2)+2*p(2)./(1+4*(w0*p(2)).^2));
lor1par_eval = @(p, w0) p(3)*(p(2)./(1+(w0*p(2)).^2)+2*p(2)./(1+4*(w0*p(2)).^2));
lor2parwou_eval = @(p,w0) p(1)*(p(2)./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2));
%3/10*(mu0/(4*pi)*gH^2*hbar/(p(1)^3))^2
%Inital guesses
p0 = [rHH0;tc0]
 
%Nonlinear fit
result = fminsearch(lor2parwou,p0)
chisqr=lor2parwou(result)
%,'MaxFunEvals',10^4,'MaxIter',10^4)
(result(1)/(3/10*((gammaH)^2*hbarSI)^2))^-1/6
(result(1)/(3/10*((gammaH)^2*hbarSI/((10^-8)^3))^2))^-1/6
plot(xdata, ydata, 'k.')
hold on
plot(xdata, lor2parwou_eval(result, xdata), 'r')
xlabel('frequency Wo')
ylabel('Relaxation rates s-1')
endfunction