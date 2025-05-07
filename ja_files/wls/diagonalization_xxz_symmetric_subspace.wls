#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory["~/Documents/chaos_meet_channels/ja_files/wls"];(*Directorio que funciona en emmy*)
Import["../../Mathematica_packages/QMB.wl"]


L = 18;
k = 5;(*subspace of k spins up*)
d = 9;(*defect site*)
Jz = 1.;


(*kernels = 30;*)
kernels = {4, 20, 36, 30, 10};


(*prueba, funciono*)
(*L = 3;
k = 2;(*subspace of k spins up*)
d = 1;(*defect site*)
Jz = 1.;*)


basis = Select[Tuples[{0,1}, L], Count[#,1] == k&];


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


SzSz=DiagonalMatrix[Total[(-1)^Differences[#]]&/@basis,TargetStructure->"Sparse"]


Sz=DiagonalMatrix[(-1)^basis[[All,d]],TargetStructure->"Sparse"]


points = Tuples[Subdivide[0.05, 5., 29], 2];


(*time1=AbsoluteTime[];*)

numkernels = Total[kernels];
kernels = Prepend[StringJoin[#1, ToString[#2]]&@@@Transpose[{{"ssh://nodo1/?","ssh://nodo2/?","ssh://nodo3/?","ssh://nodo4/?"}, kernels[[2;;]]}], kernels[[1]]];
LaunchKernels[kernels];
Print["Kernels listos: ", Length[ParallelKernels[]], "/", numkernels];

faltan = Length[points]; SetSharedVariable[faltan];
eigenvalues = 
	ParallelDo[
		Module[{\[CurlyEpsilon], Jxy, H},
		{\[CurlyEpsilon], Jxy} = i;
		t1 = AbsoluteTime[];
		H = Normal[Jxy/4*SxSxPlusSySy + Jz/4*SzSz + \[CurlyEpsilon]/2*Sz];
		Export["/home/jadeleon/Documents/chaos_meets_channels/ja_files/data/eigenenergies/xxz/L_"<>ToString[L]<>"/k_"<>ToString[k]<>"/Jz_1_omega_0_d_9_epsilon_"<>ToString[NumberForm[\[CurlyEpsilon], {Infinity, 2}]]<>"_Jxy_"<>StringPadLeft[ToString[Jxy], 2, "0"]<>".csv", Sort[Eigenvalues[H]], "CSV"];
		t2 = AbsoluteTime[];
		Print["\[CurlyEpsilon] = "<>ToString[NumberForm[\[CurlyEpsilon], {Infinity,2}]]<>", Jxy ="<>ToString[NumberForm[Jxy, {Infinity,2}]]<>": "<>ToString[NumberForm[t2 - t1, {Infinity, 0}]]<>" s. Faltan "<>ToString[faltan-=1]];
		(*Flatten[{{\[CurlyEpsilon], Jxy}, Eigenvalues[H]}]*)
		]
	, {i, points},
	DistributedContexts->Full, Method -> "FinestGrained"];
CloseKernels[];

(*time2=AbsoluteTime[];

Print["Tiempo diagonalizaciones: "<>ToString[NumberForm[time2 - time1,{Infinity,2}]]<>" s"];*)

(*Export["../data/eigenvalues_xxz_L_"<>ToString[L]<>"_d_"<>ToString[d]<>"_k_"<>ToString[k]<>"_Jz_1_omega_0.csv",Prepend[eigenvalues,{"\[CurlyEpsilon],Jxy,E_1,E_2,..."}],"CSV"];*)


Print["Termin\[EAcute]"];
