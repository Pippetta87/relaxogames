pkg load miscellaneous;
hbarSI=physical_constant('Planck constant');
gammaH=physical_constant('proton gyromag. ratio');
%Load data
%load cas13feb21;
dataribo=csvread('ribo.csv');
dataferrante=csvread('Ferrante.csv');
omegaferrante=dataferrante(:,1)*2*pi*10^6
R1alb=dataferrante(:,3)
R1H=dataferrante(:,5)
omegarib=dataribo(:,1)
R1rib=dataribo(:,2)
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

%albumine
%Create x and y values
xdata=omegaferrante*2*pi*10^6;%omegaexp
ydata=R1alb;%rexp
w0=logspace(min(xdata),max(xdata),1000);
xdatal=log(xdata);
ydatal=log(ydata);
%Inital guesses
%p0 = [25*10^9/10000, 4*10^-4];%unknown
p0_adim = [A,tc,0];
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
%semilogx(xdata, ydata, 'k.');
hold on;
%semilogx(xdata, adim(result, xdata), 'r')
%semilogx(xdata, adim(result_adim, xdata), 'r');
%semilogx(xdata, adim_lin(result_adim_lin, xdata), 'b');
xlabel('omega');
ylabel('R1alb');
hold off;
%figure;
%H2O
%Create x and y values
xdata=omegaferrante%omegaexp
ydata=R1H%rexp
w0=logspace(min(xdata),max(xdata),1000);
xdatal=log(xdata);
ydatal=log(ydata);
%Inital guesses
%p0 = [25*10^9/10000, 4*10^-4];%unknown
p0_adim = [A,tc,0];
%tc=10^-8;
%p0_adim_lin = [ydata(1)*100,tc,0,0];
%p0_adim_lin=[1.94*10^2,1.276*10^02,1.85,0];
%p0=[5,-10^-6];%exponential
%p0=[4*10^5,1];%power of x
%p0= [1,1,1];%p0k
%Nonlinear fit
options = optimset('MaxFunEvals',10000,'MaxIter',10000);
result_adim = fminsearch(@(p)adimobjective(p,xdata,ydata), p0_adim,options)
%result_adim_lin = fminsearch(@(p)adimobjective_lin(p,xdata,ydata), p0_adim_lin,options)
%semilogx(xdata, ydata, 'k.');
hold on;
%semilogx(xdata, adim(result, xdata), 'r')
%semilogx(xdata, adim(result_adim, xdata), 'r');
%semilogx(xdata, adim_lin(result_adim_lin, xdata), 'b');
xlabel('omega');
ylabel('R1H');
hold off;
%figure;
%Ribonuclease
%Create x and y values
xdata=omegarib;%omegaexp
ydata=R1rib;%rexp
w0=logspace(min(xdata),max(xdata),1000);
xdatal=log(xdata);
ydatal=log(ydata);
%Inital guesses
%p0 = [25*10^9/10000, 4*10^-4];%unknown
p0_adim = [A,tc,0];
%tc=10^-8
%p0_adim_lin = [ydata(1)*100,tc,0,0];
%p0_adim_lin=[85.7,40.6,0.5,0];
%p0=[5,-10^-6];%exponential
%p0=[4*10^5,1];%power of x
%p0= [1,1,1];%p0k
%Nonlinear fit
options = optimset('MaxFunEvals',10000,'MaxIter',10000);
result_adim = fminsearch(@(p)adimobjective(p,xdata,ydata), p0_adim,options)
%result_adim_lin = fminsearch(@(p)adimobjective_lin(p,xdata,ydata), p0_adim_lin,options)
%semilogx(xdata, ydata, 'k.');
hold on;
%semilogx(xdata, adim(result, xdata), 'r')
%semilogx(xdata, adim(result_adim, xdata), 'r');
%semilogx(xdata, adim_lin(result_adim_lin, xdata), 'b');
xlabel('omega');
ylabel('R1rib');
