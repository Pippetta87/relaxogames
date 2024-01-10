## Author: Utente <Utente@ASUS>
## Created: 2021-02-13

#omegaexp=[1.1*10^5 1.7*10^5 2.5*10^5 4.7*10^5 6*10^5 8*10^5 9*10^5 1.2*10^6 1.3*10^6 1.7*10^6 1.9*10^6 2.1*10^6 2.6*10^6 3*10^6 3.2*10^6 3.5*10^6 3.9*10^6 5.5*10^6 7.5*10^6 1.2*10^7 1.9*10^7 3.2*10^7 5.7*10^7 8.5*10^7];
#rexp=[4.47 4.45 4.40 4.38 4.16 4.1 3.8 3.6 3.4 3.1 2.8 2.5 2.2 1.9 1.7 1.4 1.2 0.8 0.52 0.32 0.18 0.13 0.10 0.08];
#freqexp2=[20 12.792 8.1769 8.1769 5.231 3.3449 2.1384 1.3682 8.74*10^-1 5.59*10^-1 3.58*10^-1 2.29*10^-1 1.46*10^-1 9.33*10^-2 5.97*10^-2 3.82*10^-2 2.46*10^-2 1.57*10^-2 9.94*10^-3];
#omegaexp2=freqexp2*2*pi;
#albumin_rexp2=[5.36*10^-1 5.62*10^-1 5.89*10^-1 6.12*10^-1 6.54*10^-1 7.36*10^-1 8.27*10^-1 9.22*10^-1 9.92*10^-1 1.0545 1.0865 1.108 1.1055 1.1227 1.1321 1.1284 1.1547 1.1394 1.1532];
#casein_d20_rexp2=[4.62*10^-1 4.83*10^-1 5.40*10^-1 5.08*10^-1 5.71*10^-1 6.16*10^-1 5.69*10^-1 6.56*10^-1 6.90*10^-1 7.84*10^-1 7.96*10^-1 9.22*10^-1 1.0105 1.1869 1.2681 1.4359 1.5096 1.6316 1.6559];
#casein_h20_rexp2=[6.72*10^-1 7.04*10^-1 7.86*10^-1 8.08*10^-1 8.69*10^-1 9.62*10^-1 1.0735 1.1681 1.3378 1.4973 1.7112 2.0102 2.3325 2.7227 3.0927 3.4644 3.8165 4.0614 4.1864];

dataribo=csvread('ribo.csv');
dataferrante=csvread('Ferrante.csv');
omegaferrante=dataferrante(:,1)*2*pi*10^6;
R1alb=dataferrante(:,3);
R1D=dataferrante(:,4);
R1H=dataferrante(:,5);
hf = figure ();
semilogx(omegaferrante,R1alb,".r;Albumin;");
hold on;
semilogx(omegaferrante,R1H,"xb;Casein in H2O;");
semilogx(omegaferrante,R1D,"^m;Casein in D2O;");
 h = legend ("location", "northeastoutside");
 print (hf, "ferranterel.png", "-dpng");
R1rib=dataribo(:,2);
omegarib=dataribo(:,1);
hf=figure;
semilogx(omegarib,R1rib,"*c;Ribonuclease;");
 print (hf, "ribo.png", "-dpng");

