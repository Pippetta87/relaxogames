pkg load miscellaneous;
hbarSI=physical_constant('Planck constant');
gammaH=physical_constant('proton gyromag. ratio');
%Load data
omegaexp=[1.1*10^5 1.7*10^5 2.5*10^5 4.7*10^5 6*10^5 8*10^5 9*10^5 1.2*10^6 1.3*10^6 1.7*10^6 1.9*10^6 2.1*10^6 2.6*10^6 3*10^6 3.2*10^6 3.5*10^6 3.9*10^6 5.5*10^6 7.5*10^6 1.2*10^7 1.9*10^7 3.2*10^7 5.7*10^7 8.5*10^7];
rexp=[4.47 4.45 4.40 4.38 4.16 4.1 3.8 3.6 3.4 3.1 2.8 2.5 2.2 1.9 1.7 1.4 1.2 0.8 0.52 0.32 0.18 0.13 0.10 0.08];
freqexp2=[20 12.792 8.1769 8.1769 5.231 3.3449 2.1384 1.3682 8.74*10^-1 5.59*10^-1 3.58*10^-1 2.29*10^-1 1.46*10^-1 9.33*10^-2 5.97*10^-2 3.82*10^-2 2.46*10^-2 1.57*10^-2 9.94*10^-3];
omegaexp2=freqexp2*2*pi;
albumin_rexp2=[5.36*10^-1 5.62*10^-1 5.89*10^-1 6.12*10^-1 6.54*10^-1 7.36*10^-1 8.27*10^-1 9.22*10^-1 9.92*10^-1 1.0545 1.0865 1.108 1.1055 1.1227 1.1321 1.1284 1.1547 1.1394 1.1532];
casein_d20_rexp2=[4.62*10^-1 4.83*10^-1 5.40*10^-1 5.08*10^-1 5.71*10^-1 6.16*10^-1 5.69*10^-1 6.56*10^-1 6.90*10^-1 7.84*10^-1 7.96*10^-1 9.22*10^-1 1.0105 1.1869 1.2681 1.4359 1.5096 1.6316 1.6559];
casein_h20_rexp2=[6.72*10^-1 7.04*10^-1 7.86*10^-1 8.08*10^-1 8.69*10^-1 9.62*10^-1 1.0735 1.1681 1.3378 1.4973 1.7112 2.0102 2.3325 2.7227 3.0927 3.4644 3.8165 4.0614 4.1864];
xdata=omegaexp;
ydata=rexp;
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

%Create x and y values
%xdata=omegaferrante*2*pi*10^6;%omegaexp
%ydata=R1alb;%rexp
w0=linspace(min(xdata),max(xdata),1000);
xdatal=log(xdata);
ydatal=log(ydata);
%Inital guesses
%p0 = [25*10^9/10000, 4*10^-4];%unknown
p0_adim = [1.1*10^7,9*10^-8,0];
%tc=2*10^-8;
%1.2384e+02   7.0207e+01   1.0219e+00  -8.5788e-16
%p0_adim_lin= [1.23*10^2,70.2,0,0];
%p0=[5,-10^-6];%exponential
%p0=[4*10^5,1];%power of x
%p0= [1,1,1];%p0k
%Nonlinear fit
options = optimset('MaxFunEvals',10000,'MaxIter',10000);
result_adim = fminsearch(@(p)adimobjective(p,xdata,ydata), p0_adim,options)
%result_adim_lin = fminsearch(@(p)adimobjective_lin(p,xdata,ydata), p0_adim_lin,options)
semilogx(xdata, ydata, 'k.');
hold on;
%semilogx(xdata, adim(result, xdata), 'r')
semilogx(xdata, adim(result_adim, xdata), 'r');
%semilogx(xdata, adim_lin(result_adim_lin, xdata), 'b');
xlabel('omega');
ylabel('R1alb');
hold off;
