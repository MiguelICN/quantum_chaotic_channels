#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Anotaciones*)


(* ::Section:: *)
(*Definiciones para este script*)


(* ::Section:: *)
(*Par\[AAcute]metros*)


(* ::Text:: *)
(*Parametros del sistema*)


L = 7;(* espines *)
d = 3;(* defecto *)
Jz = 1;
\[Omega] = 0;


(* ::Text:: *)
(*Tiempos*)


{t0, tf, deltat}={0., 49.5, 0.5};(* t inicial, t final, n\[UAcute]mero de puntos *)


(* ::Text:: *)
(*Numero de kernels a utilizar en la paralelizacion*)


(*kernels = 24;*)


(* ::Section:: *)
(*Ac\[AAcute] comienza el c\[AAcute]lculo*)


(* ::Text:: *)
(*Fijar directorio, cargar paquete QMB y encender todos los kernels a utilizar*)


SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files"];
Get["../Mathematica_packages/QMB.wl"];


(* ::Text:: *)
(*Inicializar variables*)


tiempos = Table[t, {t, t0, tf, deltat}];

points = Tuples[Subdivide[0.05, 5., 29], 2][[125;;]];
lenpoints = Length[points];

(* Pauli strings *)
pauliIndices = Tuples[Range[0,3],L-1];
pauliStrings = Pauli /@ pauliIndices;

(* Constantes que aparecen dentro de la suma del promedio de Haar *)
constants = (1/2)^2*(1/12)^(L-1)*Power[3., (L-1) - Count[#, x_/; x!=0]] & /@ pauliIndices;

(* Lista de tiempos por cada iteraci\[OAcute]n del Do *)
times={};

dimH=2^L;
halfDimH=dimH/2;
halfDimHPlusOne=dimH/2+1;


LaunchKernels[25];
Print["N\[UAcute]mero de kernels listos: ", Length[ParallelKernels[]]];


Do[
startTime = AbsoluteTime[];

(* Asignar variables \[CurlyEpsilon] y Jxy *)
{\[CurlyEpsilon], Jxy} = points[[i]];

(* Crear Hamiltoniano y diagonalizar *)
H = XXZOpenHamiltonian[Jxy, Jz, \[Omega], \[CurlyEpsilon], L, d];
{eigenvalues, eigenvectors} = Chop[Eigensystem[H]];

(* Calcular la matriz de cambio de base de la de eigenenerg\[IAcute]as a la computacional *)
P = Transpose[eigenvectors];

(* Operador de evoluci\[OAcute]n *)
ClearAll[U]; U[t_] := P . DiagonalMatrix[Exp[-I*eigenvalues*t]] . ConjugateTranspose[P];

(* Crear archivo para guardar valor esperado de Haar de la pureza de Choi *)
(*HaarAvgChoiPurityFilename = "data/haar_avg_choi_purity/xxz/L_"<>ToString[L]<>"/d_"<>ToString[d]<>"_omega_"<>ToString[\[Omega]]<>"_Jz_1_epsilon_"<>ToString[NumberForm[\[CurlyEpsilon],{Infinity,2}]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>".csv";
dataFile = OpenAppend[HaarAvgChoiPurityFilename];
WriteString[dataFile, "'t', 'valor esperado de Haar de la pureza de Choi'\n"];
Close[dataFile];*)

(* Calcular valor esperado de Haar de la pureza de Choi para cada t *)
DistributeDefinitions[U, Pauli, L];
choipurity = {}; SetSharedVariable[choipurity];
ParallelDo[
	t1 = AbsoluteTime[];
	(* Asignar el operador de evoluci\[OAcute]n a una variable (para no repetir en cada t\[EAcute]rmino de la suma del valor esperado de Haar) *)
	u = U[t]; uDg = ConjugateTranspose[u];
	{{u11, u12}, {u21, u22}} = {{u[[1;;halfDimH, 1;;halfDimH]], u[[1;;halfDimH, halfDimHPlusOne;;dimH]]}, {u[[halfDimHPlusOne;;dimH, 1;;halfDimH]], u[[halfDimHPlusOne;;dimH, halfDimHPlusOne;;dimH]]}};
	{{uDg11, uDg12}, {uDg21, uDg22}} = {{uDg[[1;;halfDimH, 1;;halfDimH]], uDg[[1;;halfDimH, halfDimHPlusOne;;dimH]]}, {uDg[[halfDimHPlusOne;;dimH, 1;;halfDimH]], uDg[[halfDimHPlusOne;;dimH, halfDimHPlusOne;;dimH]]}};
	(*DistributeDefinitions[u11, u12, u21, u22, uDg11, uDg12, uDg21, uDg22, Pauli];*)
	haarAvgChoiPurity = 
	Total[
		Table[
		With[{c = constants[[i]]},
			(* Calcular \Tr_S( U(t).\sigma .U^\dagger(t) ) *)
			A = u11 . # . uDg11 + u12 . # . uDg21 + u21 . # . uDg12 + u22 . # . uDg22 &[Pauli[IntegerDigits[i-1, 4, L - 1]]];
			(* Calcular la pureza de A y multiplicarlo por su respectiva constante *)
			Chop[c * Total[Abs[A]^2, 2]]
		]
		, {i, 4^(L-1)}]
	];
	t2 = AbsoluteTime[]; Print["(\[CurlyEpsilon]="<>ToString[NumberForm[\[CurlyEpsilon], {Infinity, 2}]]<>", Jxy="<>ToString[NumberForm[Jxy, {Infinity, 2}]]<>"). t = "<>ToString[NumberForm[t, {Infinity, 2}]]<>": "<>ToString[t2 - t1]<>" s"];
	AppendTo[choipurity, {t, haarAvgChoiPurity}];
	(* Agregar {t, valor esperado de Haar de la pureza de Choi} al archivo CSV *)
	(*dataFile = OpenAppend[HaarAvgChoiPurityFilename];
	WriteString[dataFile, StringJoin[Riffle[ToString/@{t, haarAvgChoiPurity}, ","]]<>"\n"];
	Close[dataFile];*)
,{t, tiempos}];
(* Exportar datos *)
filename = "data/haar_avg_choi_purity/xxz/L_"<>ToString[L]<>"/d_"<>ToString[d]<>"_omega_"<>ToString[\[Omega]]<>"_Jz_1_epsilon_"<>ToString[NumberForm[\[CurlyEpsilon],{Infinity,2}]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>".csv";
Export[filename, Prepend[SortBy[choipurity, First], {"t", "valor esperado de Haar de la pureza de Choi"}], "CSV"];

endTime = AbsoluteTime[];
time = endTime-startTime;
remainingIterations = (lenpoints - i);
AppendTo[times, time];
Print[DateString[], ". ", "(\[CurlyEpsilon]="<>ToString[NumberForm[\[CurlyEpsilon], {Infinity, 2}]]<>", Jxy="<>ToString[NumberForm[Jxy, {Infinity, 2}]]<>") done. ", ToString[i], "/", lenpoints,": ", ToString[NumberForm[time,{Infinity,0}]], " s. Terminar\[AAcute]: ", DateString[DatePlus[Now,Quantity[remainingIterations*Mean[times]/3600,"Hours"]]]];

, {i, lenpoints}]; 
