#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory["~/Documents/chaos_meet_channels/ja_files/wls"];(*Directorio que funciona en emmy*)
Import["../../Mathematica_packages/QMB.wl"]


L = 18;
k = 7;(*subspace of k spins up*)
d = 9;(*defect site*)
Jz = 1.;


(*Table[{k,L!/(k!(L-k)!)},{k,L}]*)


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


(*LaunchKernels[3];
Print["Kernels encedidos: ", Length[ParallelKernels[]]];
AbsoluteTiming[
SxSxPlusSySy=Module[
	{basisSubspace=Tuples[FromDigits[#,2]+1&/@basis,2],
	dimsubspace=Length[basis],
	SxSxPlusSySy=With[{NNindices = RotateRight[Join[{1,1},ConstantArray[0,L-2]],#]&/@Range[0,L-2]},
	N[Total[Join[Pauli/@NNindices,Pauli/@(2NNindices)]]]],indicesones,intersection
	},
	indicesones=SxSxPlusSySy["ExplicitPositions"];
	intersection=Intersection[basisSubspace,indicesones];
	DistributeDefinitions[basisSubspace, SxSxPlusSySy, intersection, dimsubspace];
	ArrayReshape[SparseArray[
		ParallelMap[
			Flatten[Position[basisSubspace,#]]->SxSxPlusSySy[[#[[1]],#[[2]]]]&, intersection,
			Method -> "CoarsestGrained"
		]
	,dimsubspace^2],{dimsubspace,dimsubspace}]
];
]
CloseKernels[];*)


time2=AbsoluteTime[];
Print["Construcci\[OAcute]n de SxSxPlusSySy tom\[OAcute] "<>ToString[NumberForm[time2-time1,{Infinity,2}]]<>" s"];


SzSz=DiagonalMatrix[Total[(-1)^Differences[#]]&/@basis,TargetStructure->"Sparse"]


Sz=DiagonalMatrix[(-1)^basis[[All,d]],TargetStructure->"Sparse"]


points = Tuples[Subdivide[0.05, 5., 29], 2];


time1=AbsoluteTime[];

LaunchKernels[30];
eigenvalues = 
	ParallelTable[
		Module[{\[CurlyEpsilon], Jxy, H},
		{\[CurlyEpsilon], Jxy} = i;
		H = Normal[Jxy/4*SxSxPlusSySy + Jz/4*SzSz + \[CurlyEpsilon]/2*Sz];
		Flatten[{{\[CurlyEpsilon], Jxy}, Eigenvalues[H]}]
		]
	, {i, points},
	DistributedContexts->Full, Method->"CoarsestGrained"];
CloseKernels[];

time2=AbsoluteTime[];

Print["Tiempo diagonalizaciones: "<>ToString[NumberForm[time2 - time1,{Infinity,2}]]<>" s"];

Export["../data/eigenvalues_xxz_L_"<>ToString[L]<>"_d_"<>ToString[d]<>"_k_"<>ToString[k]<>"_Jz_1_omega_0.csv",Prepend[eigenvalues,{"\[CurlyEpsilon],Jxy,E_1,E_2,..."}],"CSV"];
