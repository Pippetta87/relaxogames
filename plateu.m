%function r = g (x)
%  r = [ sumsq(x)-10;
%        x(2)*x(3)-5*x(4)*x(5);
%        x(1)^3+x(2)^3+1 ];
%endfunction
%function obj = phi (x)
%  obj = exp (prod (x)) - 0.5*(x(1)^3+x(2)^3+1)^2;
%endfunction
%x0 = [-1.8; 1.7; 1.9; -0.8; -0.8];
%[x, obj, info, iter, nf, lambda] = sqp (x0, @phi, @g, [])
%x =
%  -1.71714
%   1.59571
%   1.82725
%  -0.76364
%  -0.76364
%obj = 0.053950
%info = 101
%iter = 8
%nf = 10
%lambda =
%-0.0401627
%   0.0379578
%  -0.0052227
clear all
pkg load miscellaneous;
hbarSI=physical_constant('Planck constant');
gammaH=physical_constant('proton gyromag. ratio');
dataribo=csvread('ribo.csv');
dataferrante=csvread('Ferrante.csv');
omegaferrante=dataferrante(:,1)*2*pi*10^6
R1alb=dataferrante(:,3)
R1D=dataferrante(:,4)
R1H=dataferrante(:,5)
omegarib=dataribo(:,1)
R1rib=dataribo(:,2)
%Define the exponential model
%cgs: gauss=10^4T mu0=1 J=10^7erg
Rcgs = @(param,x) 3/10*(gammaH*10^-4)^2*(10^7*hbarSI)/(param(1)^3*(10^-8)^3)^2*param(2).*(1/(1+(x*param(2)).^2)+2./(1+4*(x*param(2)).^2));
cgs_objective=@(param,x,y) sum((y-adim(param,x)).^2./y);
adim = @(param,x) param(1)*(param(2)./(1+x.^2*param(2)^2)+4*param(2)./(1+4*x.^2*param(2)^2))+param(3);
%adimobjective = @(param,x,y) sum((y-adim(param,x)).^2./y);
adim_lin = @(param,x) param(1)*(param(2)./(1+x.^2*param(2)^2)+4*param(2)./(1+4*x.^2*param(2)^2))+param(3)+param(4)*x;
adimobjective_lin = @(param,x,y) sum((y-adim_lin(param,x)).^2./y);
global xdata=omegaferrante
global ydata=R1D
U=find(abs(ydata-max(ydata))<0.2)
L=find(abs(ydata-min(ydata))<0.2)
   chiquadro=@(p) sum((adim(p,xdata)).^2./ydata);
%X=fminsearch(chiquadro,[round(length(obscut)/2);1;5000;100])
global datalength=length(xdata)
function zeroing=P2P(p)
global datalength
global ydata
global xdata
adim = p(1)*(p(2)./(1+xdata.^2*p(2)^2)+4*p(2)./(1+4*xdata.^2*p(2)^2))+p(3);
zeroing([1:datalength])=(ydata([1:datalength])-adim([1:datalength])).^2;
zeroing(datalength+1)=max(10^5-p(1),0);
zeroing(datalength+2)=max(-10^9+p(1),0);
zeroing(datalength+3)=max(-10^-6+p(2),0);
zeroing(datalength+4)=max(10^-9-p(2),0);
zeroing(datalength+5)=max(-p(2),0)
endfunction
fh=figure;
plot(xdata,ydata,'.')
hold on
p0_adim = [0.5*10^7,2*10^-7,0];
options.TolFun=10^-40;options.TolX=10^-40;
%options = optimset('MaxFunEvals',10000,'MaxIter',10000);
%result_adim = fminsearch(@(p)adimobjective(p,xdata(U),ydata(U)), p0_adim,options)
%result_adim = fminsearch(@(p)adimobjective(p,xdata(L),ydata(L)), p0_adim,options)
[relaxpar,info]=fsolve("P2P", [10^6 10^-8 1]',options);%vector of param
semilogx(xdata, adim(relaxpar, xdata), 'r');
legend("$R_1^{NS}$ misurato", "$R_1^{NS}(\\tau_c)$",'Interpreter','latex')
saveas(fh,"R1D.png");
