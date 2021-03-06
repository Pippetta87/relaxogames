## Copyright (C) 2021 Utente
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {} {@var{retval} =} caseinate (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: Utente <Utente@ASUS>
## Created: 2021-02-13

#omegaexp=[1.1*10^5 1.7*10^5 2.5*10^5 4.7*10^5 6*10^5 8*10^5 9*10^5 1.2*10^6 1.3*10^6 1.7*10^6 1.9*10^6 2.1*10^6 2.6*10^6 3*10^6 3.2*10^6 3.5*10^6 3.9*10^6 5.5*10^6 7.5*10^6 1.2*10^7 1.9*10^7 3.2*10^7 5.7*10^7 8.5*10^7];
#rexp=[4.47 4.45 4.40 4.38 4.16 4.1 3.8 3.6 3.4 3.1 2.8 2.5 2.2 1.9 1.7 1.4 1.2 0.8 0.52 0.32 0.18 0.13 0.10 0.08];
yah=[0.54 0.57 0.62 0.66 0.70 0.78 0.89 0.98 1.07 1.09 1.11 1.11 1.12 1.12 1.12 1.12 1.12 1.13 1.13];
xah=[1.26*10^8 8.034*10^7 5.135*10^7 3.286*10^7 2.10*10^7 1.34*10^7 8.59*10^6 5.49*10^6 3.51*10^6 2.25*10^6 1.44*10^6 9.17*10^5 5.87*10^5 3.75*10^5 2.4*10^5 1.54*10^5 9.86*10^4 6.24*10^4];
ych=[0.22 0.24 0.26 0.30 0.37 0.49 0.64 1.00 1.32 1.75 2.30 2.75 3.19 3.56 3.88 4.07 4.19 4.25 4.28];
xch=[1.26*10^8 8.034*10^7 5.135*10^7 3.286*10^7 2.10*10^7 1.34*10^7 8.59*10^6 5.49*10^6 3.51*10^6 2.25*10^6 1.44*10^6 9.17*10^5 5.87*10^5 3.75*10^5 2.4*10^5 1.54*10^5 9.86*10^4 6.24*10^4 4.05*10^4];
ycd=[0.48 0.50 0.53 0.53 0.57 0.59 0.59 0.67 0.71 0.80 0.96 1.14 1.37 1.57 1.63 1.65 1.66 1.68 1.70];
xcd=[1.26*10^8 8.034*10^7 5.135*10^7 3.286*10^7 2.10*10^7 1.34*10^7 8.59*10^6 5.49*10^6 3.51*10^6 2.25*10^6 1.44*10^6 9.17*10^5 5.87*10^5 3.75*10^5 2.4*10^5 1.54*10^5 9.86*10^4 6.24*10^4];