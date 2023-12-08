hbarSI=physical_constant('Planck constant');
gammaH=physical_constant('proton gyromag. ratio');
%Load data
%load cas13feb21;
%Create x and y values
xdata=omegaexp2;%omegaexp
ydata=casein_h20_rexp2;%rexp
w0=linspace(min(xdata),max(xdata),1000);
xdatal=log(omegaexp);
ydatal=log(rexp);
%Define the exponential model
%cgs: gauss=10^4T mu0=1 J=10^7erg

my_fun = @(p) sum((ydata-3/10*((gammaH*10^-4)^2*(10^7*hbarSI)/(p(1)^3*(10^-8)^3))^2*(p(2)./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2./ydata);
my_fun_two = @(p) sum((ydata-3/10*gammaH^2*hbarSI/p(1)^6*p(2)./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2)).^2./ydata);
my_funl = @(p) sum((ydatal-log(3/10*(gammaH*10^-4)^2*10^7*hbarSI)^2+6*log(p(1)*(10^-8)^6)-log(p(2))+log(1./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2);
my_funk = @(p) sum((ydata-p(3)*(p(2)./(1+(xdata*p(2)).^2)+2*p(2)./(1+4*(xdata*p(2)).^2))).^2/ydata);

%Function to evaluate the model. In this case we need to provide the error
%funtion: error sum of squares
fun_eval = @(p, w0) 3/10*((gammaH*10^-4)^2*hbarSI*10^7/(p(1)^3*(10^-8)^3))^2*(p(2)./(1+(w0*p(2)).^2)+2*p(2)./(1+4*(w0*p(2)).^2));
fun_eval_two = @(p, w0) 3/10*gammaH^2*hbarSI/p(1)^6*(p(2)./(1+(w0*p(2)).^2))+2*p(2)./(1+4*(w0*p(2)).^2);
fun_evall = @(p, w0) log(3/10)+2*log((10^-4*gammaH)^2*hbarSI*10^7)-6*log(p(1)*(10^-8)^6)+log(p(2))-log(1./(1+(w0*p(2)).^2)+2*p(2)./(1+4*(w0*p(2)).^2));
fun_evalk = @(p, w0) p(3)*(p(2)./(1+(w0*p(2)).^2)+2*p(2)./(1+4*(w0*p(2)).^2));
%3/10*(mu0/(4*pi)*gH^2*hbar/(p(1)^3))^2

adim = @(param,x) param(1)*(param(2)./(1+x.^2*param(2)^2)+4*param(2)./(1+4*x.^2*param(2)^2))+param(3);
adim_lin = @(param,x) param(1)*(param(2)./(1+x.^2*param(2)^2)+4*param(2)./(1+4*x.^2*param(2)^2))+param(3)+param(4)*x;
exponential = @(param,x) param(1)*exp(param(2)*x);
potenza=@(param,x) param(1)./x.^param(2);
adimobjective = @(param,x,y) sum((y-adim(param,x)).^2./y);
adimobjective_lin = @(param,x,y) sum((y-adim_lin(param,x)).^2./y);

%Inital guesses
%p0 = [25*10^9/10000, 4*10^-4];%unknown
p0_adim = [0.045,1/0.225,0];
p0_adim_lin = [0.045,1/0.225,0,0];
%p0=[5,-10^-6];%exponential
%p0=[4*10^5,1];%power of x
%p0= [1,1,1];%p0k
%Nonlinear fit
options = optimset('MaxFunEvals',10000,'MaxIter',10000)
result_adim = fminsearch(@(p)adimobjective(p,xdata,ydata), p0_adim,options)
result_adim_lin = fminsearch(@(p)adimobjective_lin(p,xdata,ydata), p0_adim_lin,options)
plot(xdata, ydata, 'k.')
hold on
%semilogx(xdata, adim(result, xdata), 'r')
semilogx(xdata, adim(result_adim, xdata), 'r')
semilogx(xdata, adim_lin(result_adim_lin, xdata), 'b')
xlabel('omega')
ylabel('R')
