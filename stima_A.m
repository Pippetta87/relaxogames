angles=[38.59,38.40,7.62,7.21]
angles_rad=angles*pi/180
tc=[3.0*10^-8, 6.16*10^-7, 6.25*10^-7]
M=(tc*10^9-0.1674)/0.0005998
load pkg miscellaneous
pkg load miscellaneous
Na=physical_constant("Avogadro constant")
M_mg=M/Na*1000
ang_coeff=tan(angles_rad)*graph_conv
A=ang_coeff.*M_mg
