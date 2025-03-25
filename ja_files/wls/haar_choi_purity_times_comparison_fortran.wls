#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Definiciones para este script*)


(* Funci\[OAcute]n para tiempos espaciados en escala log *)
LogSpace[start_,end_,n_] := Exp@Subdivide[Log[start],Log[end],n-1];


(* ::Section:: *)
(*Par\[AAcute]metros (s\[OAcute]lo modificar esta secci\[OAcute]n)*)


(* ::Text:: *)
(*Parametros del sistema*)


L = 10;(* espines *)
hx = 1.;

{hz0, hzf}={0.05, 2.95};(* se considera un grid de n^2 puntos equiespaciados *)
n = 10;(* n\[UAcute]mero de puntos en el mapa *)


(* ::Text:: *)
(*Tiempos*)


{t0, tf, tn}={0.1, 1000., 100};(* t inicial, t final, n\[UAcute]mero de puntos *)


(* ::Text:: *)
(*Numero de kernels a utilizar en la paralelizacion*)


kernels= 1;


(* ::Section:: *)
(*Ac\[AAcute] comienza el c\[AAcute]lculo*)


(* ::Text:: *)
(*Fijar directorio, cargar paquete QMB y encender todos los kernels a utilizar*)


SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files"];
Get["../Mathematica_packages/QMB.wl"];
(*LaunchKernels[kernels];*)
Print["N\[UAcute]mero de kernels listos: ", Length[ParallelKernels[]]];


(* ::Text:: *)
(*Inicializar variables*)


tiempos = LogSpace[t0, tf, tn];

points = Tuples[{#,#}]&[Subdivide[hz0, hzf, n - 1]];
lenpoints = Length[points];

(* Pauli strings *)
pauliIndices = Tuples[Range[0,3],L-1][[2;;2]];
pauliStrings = Pauli /@ pauliIndices;

(* Constantes que aparecen dentro de la suma del promedio de Haar *)
constants = (1/2)^2*(1/12)^(L-1)*Power[3., (L-1) - Count[#, x_/; x!=0]] & /@ pauliIndices;

(* Lista de tiempos por cada iteraci\[OAcute]n del Do *)
times={};

dimH=2^L;
halfDimH=dimH/2;
halfDimHPlusOne=dimH/2+1;


Do[
startTime = AbsoluteTime[];

(* Asignar variables hz y J *)
{hz, J} = points[[hzAndJ]];

(* Crear Hamiltoniano y diagonalizar *)
H = IsingNNOpenHamiltonian[hx, hz, J, L];
{eigenvalues, eigenvectors} = Chop[Eigensystem[H]];

(* Calcular la matriz de cambio de base de la de eigenenerg\[IAcute]as a la computacional *)
P = Transpose[eigenvectors];

(* Operador de evoluci\[OAcute]n *)
ClearAll[U]; U[t_] := P . DiagonalMatrix[Exp[-I*eigenvalues*t/J]] . ConjugateTranspose[P];

(* Asignar el operador de evoluci\[OAcute]n a una variable (para no repetir en cada t\[EAcute]rmino de la suma del valor esperado de Haar) *)
With[{t=1.3},
	u = U[t]; uDg = ConjugateTranspose[u];
	{{u11, u12}, {u21, u22}} = {{u[[1;;halfDimH, 1;;halfDimH]], u[[1;;halfDimH, halfDimHPlusOne;;dimH]]}, {u[[halfDimHPlusOne;;dimH, 1;;halfDimH]], u[[halfDimHPlusOne;;dimH, halfDimHPlusOne;;dimH]]}};
	{{uDg11, uDg12}, {uDg21, uDg22}} = {{uDg[[1;;halfDimH, 1;;halfDimH]], uDg[[1;;halfDimH, halfDimHPlusOne;;dimH]]}, {uDg[[halfDimHPlusOne;;dimH, 1;;halfDimH]], uDg[[halfDimHPlusOne;;dimH, halfDimHPlusOne;;dimH]]}};
]

,{hzAndJ, 1}]; 


Export["pauliString.txt", pauliStrings[[1]], "Table"]


(* Convertir la matriz en pares {Re, Im} y exportar *)
Export["pauliString.txt", 
  Flatten[Map[Function[z, {Re[z], Im[z]}], pauliStrings[[1]], {2}],1], 
  "Table"]



AbsoluteTiming[
A = u11 . # . uDg11 + u12 . # . uDg21 + u21 . # . uDg12 + u22 . # . uDg22 &[ pauliStrings[[1]] ];
Chop[Total[Abs[A]^2, 2]]]
