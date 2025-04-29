#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Notas*)


(* ::Text:: *)
(*- cambiar l\[IAcute]nea del directorio, ya *)
(*- cambiar k, ya*)
(*- cambiar el numero de kernels en la paralelizacion, ya*)


(* ::Section:: *)
(*Codigo*)


SetDirectory["~/Documents/chaos_meets_channels/ja_files/wls"];(*jungkook, rm y robot*)
(*SetDirectory["~/Documents/chaos_meet_channels/ja_files/wls"];(*emmy*)*)
Import["../../Mathematica_packages/QMB.wl"]


L = 18;
k = 5;(*subspace of k spins up*)


(*Table[{k,L!/(k!(L-k)!)},{k,L}]*)


(*(*prueba, funciono*)
L = 3;
k = 2;(*subspace of k spins up*)*)


basis = Select[Tuples[{0,1}, L], Count[#, 1] == k &];


time1=AbsoluteTime[];


SxSxPlusSySy=Module[
	{SxSxPlusSySy=
	With[{NNindices = RotateRight[Join[{1,1},ConstantArray[0,L-2]],#]&/@Range[0,L-2]},
		N[Total[Join[Pauli/@NNindices,Pauli/@(2NNindices)]]]
	],
	l = FromDigits[#,2] &/@ basis + 1
	},
	SxSxPlusSySy[[l,l]]
];


time2=AbsoluteTime[];
Print["Construcci\[OAcute]n de SxSxPlusSySy tom\[OAcute] "<>ToString[NumberForm[time2-time1,{Infinity,2}]]<>" s"];


SzSz=DiagonalMatrix[Total[(-1)^Differences[#]]&/@basis,TargetStructure->"Sparse"];


points = Flatten[Table[{i, j} ,{i, Subdivide[0.05, 5, 29]}, {j, 30}], 1];


time1=AbsoluteTime[];

LaunchKernels[30];(*emmy, jungkook, r2d2*)
(*LaunchKernels[{"ssh://nodo1/?4","ssh://nodo2/?4","ssh://nodo3/?20","ssh://nodo4/?10"}](*robot*)*)
Print["Kernels: ", Length[ParallelKernels[]], "/15"];

eigenvalues = 
	ParallelTable[
		Module[{h, j, hlist, Sz, H, rand},
		{h, j} = i;(* j indexa la realizaci\[OAcute]n de h *)
		rand = SeedRandom[j*30284];Print["Random number: ", j*30284];
		hlist = RandomReal[{-h, h}, L];
		Print["hlist =  ", NumberForm[hlist, {Infinity, 2}]];
		Sz = Total[hlist[[#]]*DiagonalMatrix[(-1)^basis[[All, #]], TargetStructure->"Sparse"]&/@ Range[L]];
		H = Normal[1/4*(SxSxPlusSySy + SzSz) + 1/2*Sz];
		Export["/home/jadeleon/Documents/chaos_meets_channels/ja_files/data/eigenenergies/heisenberg/L_"<>ToString[L]<>"/k_"<>ToString[k]<>"/h_"<>ToString[NumberForm[h, {Infinity, 2}]]<>"_"<>StringPadLeft[ToString[j], 2, "0"]<>".csv", Flatten[{h, Sort[Eigenvalues[H]]}], "CSV"]
		(*Export["../data/eigenenergies/heisenberg/L_"<>ToString[L]<>"/k_"<>ToString[k]<>"/h_"<>ToString[NumberForm[h, {Infinity, 2}]]<>"_"<>StringPadLeft[ToString[j], 2, "0"]<>".csv", Flatten[{h, Sort[Eigenvalues[H]]}], "CSV"]*)
		]
	, {i, points},
	DistributedContexts->Full, Method->"CoarsestGrained"];
	
CloseKernels[];

time2=AbsoluteTime[];

Print["Tiempo diagonalizaciones: "<>ToString[NumberForm[time2 - time1,{Infinity,2}]]<>" s"];


(*eigenenergies = Table[
	Module[{h, j},
		{h, j} = i;(* j indexa la realizaci\[OAcute]n de h *)
		Flatten[Import["/home/jadeleon/Documents/chaos_meets_channels/ja_files/data/eigenenergies/heisenberg/L_"<>ToString[L]<>"/k_"<>ToString[k]<>"/h_"<>ToString[NumberForm[h, {Infinity, 2}]]<>"_"<>StringPadLeft[ToString[j], 2, "0"]<>".csv", "CSV"]]
		(*Flatten[Import["../data/eigenenergies/heisenberg/L_"<>ToString[L]<>"/k_"<>ToString[k]<>"/h_"<>ToString[NumberForm[h, {Infinity, 2}]]<>"_"<>StringPadLeft[ToString[j], 2, "0"]<>".csv", "CSV"]]*)
	]
, {i, points}];

Export["../data/eigenenergies/heisenberg/L_"<>ToString[L]<>"_k_"<>ToString[k]<>".csv", Prepend[eigenenergies, "h, E_1, E_2, ... "], "CSV"];*)
