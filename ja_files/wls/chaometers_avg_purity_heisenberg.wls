#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Set directory to this notebook's directory *)
SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files"];
Get["../Mathematica_packages/QMB.wl"];
Print["Paquetes cargados. ", DateString[],"."];


(* Par\[AAcute]metros fijos L *)
L = 7; 
{t0, tf, tn}={0., 50., 0.5};(* t inicial, t final, n\[UAcute]mero de puntos *)
kernels = 25;


(* Lista de tiempos por cada iteraci\[OAcute]n del Do *)
times = {};

points = Flatten[Table[{i, j} , {j, 30}, {i, Subdivide[0.05, 5, 29]}], 1];
lenpoints = Length[points];


(* Encender todos los kernels *)
LaunchKernels[kernels];
Print["N\[UAcute]mero de kernels disponibles: ", Length[ParallelKernels[]]];
(* ------------------------- *)

Do[
startTime = AbsoluteTime[];

{h, j} = points[[i]];(* j indexa la realizaci\[OAcute]n de h *)
hlist = RandomReal[{-h, h}, L];

(* Crear Hamiltoniano y diagonalizar *)
H = HeisenbergXXXwNoise[hlist, L];
{eigenvalues, eigenvectors} = Chop[Eigensystem[H]];
(*-----------------------------------------------------*)

(* Calcular la curva (t, pureza del caometro) *)
\[Psi]0 = RandomChainProductState[L];
purityovertime =
	ParallelTable[
		\[Psi] = Chop[StateEvolution[t, \[Psi]0, eigenvalues, eigenvectors]];
		\[Rho]A = MatrixPartialTrace[Dyad[\[Psi]], 2, {2, 2^(L-1)}];
		{t, Chop[Purity[\[Rho]A]]}
	, {t, t0, tf, tn}, DistributedContexts -> Full, Method -> "CoarsestGrained"];
(*-----------------------------------------------------*)
			
(* Exportar datos de la curva de pureza *)
filename = "data/chaometer_purity/heisenberg/L_"<>ToString[L]<>"/h_"<>ToString[NumberForm[h, {Infinity, 2}]]<>"_"<>StringPadLeft[ToString[j], 3, "0"]<>".csv";
Export[filename, Prepend[purityovertime, {"t", "pureza espin en el extremo izquierdo"}], "Table"];
(*-----------------------------------------------------*)

endTime = AbsoluteTime[];
time = endTime - startTime;
remainingIterations = (lenpoints - i);
AppendTo[times, time];
Print[DateString[], ". ",ToString[i], "/", lenpoints,": ", ToString[NumberForm[time,{Infinity,0}]], " s. Terminar\[AAcute]: ", DateString[DatePlus[Now,Quantity[remainingIterations*Mean[times]/3600,"Hours"]]]];

, {i, lenpoints}];

CloseKernels[];


Print["Termin\[EAcute], ", DateString[],"."];
