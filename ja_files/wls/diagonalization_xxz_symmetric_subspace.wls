#!/usr/bin/env wolframscript
(* ::Package:: *)

Import["../../Mathematica_packages/QMB.wl"]


(*L = 18;
k = 7;(*subspace of k spins up*)
d = 9;(*defect site*)
Jz = 1.;*)


L = 3;
k = 2;(*subspace of k spins up*)
d = 1;(*defect site*)
Jz = 1.;


LaunchKernels[30];
Print["Kernels encedidos: ", Length[ParallelKernels[]]];


basis = Select[Tuples[{0,1}, L], Count[#,1] == k&];


time1=AbsoluteTime[];


SxSxPlusSySy=Module[
	{basisSubspace=Tuples[FromDigits[#,2]+1&/@basis,2],
	dimsubspace=Length[basis],
	SxSxPlusSySy=With[{NNindices = RotateRight[Join[{1,1},ConstantArray[0,L-2]],#]&/@Range[0,L-2]},
	N[Total[Join[Pauli/@NNindices,Pauli/@(2NNindices)]]]],indicesones,intersection
	},
	indicesones=SxSxPlusSySy["ExplicitPositions"];
	intersection=Intersection[basisSubspace,indicesones];
	ArrayReshape[SparseArray[ParallelMap[Flatten[Position[basisSubspace,#]]->SxSxPlusSySy[[#[[1]],#[[2]]]]&,intersection],dimsubspace^2],{dimsubspace,dimsubspace}]
]


time2=AbsoluteTime[];
Print["Construcci\[OAcute]n de SxSxPlusSySy tom\[OAcute] "<>ToString[NumberForm[time2-time1,{Infinity,2}]]<>" s"];


SzSz=DiagonalMatrix[Total[(-1)^Differences[#]]&/@basis,TargetStructure->"Sparse"]


Sz=DiagonalMatrix[(-1)^basis[[All,d]],TargetStructure->"Sparse"]


points = Tuples[Subdivide[0.05, 5., 29],2];


time1=AbsoluteTime[];

eigenvalues = 
	ParallelTable[
		Module[{\[CurlyEpsilon], Jxy, H},
		{\[CurlyEpsilon], Jxy} = i;
		H = Normal[Jxy*SxSxPlusSySy + Jz*SzSz + \[CurlyEpsilon]*Sz];
		Flatten[{{\[CurlyEpsilon], Jxy}, Eigenvalues[H]}]
		]
	, {i, points},
	DistributedContexts->Full, Method->"CoarsestGrained"];
CloseKernels[];

time2=AbsoluteTime[];

Print["Tiempo diagonalizaciones: "<>ToString[NumberForm[time2 - time1,{Infinity,2}]]<>" s"];

Export["data_eigenvalues_xxz_L_"<>ToString[L]<>"_Jz_1_omega_0.csv",Prepend[eigenvalues,{"\[CurlyEpsilon],Jxy,E_1,E_2,..."}],"CSV"];


(* ::Subsection::Closed:: *)
(*Nueva idea*)


(*RepeatedTiming[ConstantArray[0.,{2^#,2^#}];]&/@Range[3,10]*)


(*base2=DeleteDuplicates[If[PalindromeQ[#],{#},Sort[{#,Reverse[#]}]]&/@Tuples[{0,1},L]]*)


(*Map[Identity[3]-ConstantArray[#, 3]&,base2,{2}]*)


(*(*tags*)
L = 3;

base = PositiveParitySubspaceBasis[L];
tagsbase = Tag /@ base;

AbsoluteTiming[
(*Esto seguro se puede optimizar, pero queda para despu\[EAcute]s*)
Hx = Module[{m ,dimSubspace = Length[tagsbase], q, l, tmp1},
m = SparseArray[Prepend[
	Table[
		(* Encontrar estados \sum_j \dyad{0}{1}_j \ket{k} *)
		tmp1 = Select[Join[ConstantArray[j,L]-IdentityMatrix[L],If[PalindromeQ[j],{},ConstantArray[Reverse[j],L]-IdentityMatrix[L]]], And@@NonNegative[#]&];
		q = ReplaceAll[PalindromeQ[j],{True->0,False->1}];
		(* Encontrar la posici\[OAcute]n en la base y asignar un valor *)
		Total[
			Table[
				l=ReplaceAll[PalindromeQ[k], {True->0, False->1}];
				SparseArray[Flatten[Position[tagsbase, If[MemberQ[tagsbase, Tag[k]], Tag[k], Tag[Reverse[k]]]]]-> 2.^(-(q+l)/2), dimSubspace]
		, {k, tmp1}]
		]
	,{j,base[[2;;]]}]
	,SparseArray[{},dimSubspace]
	]
	];
Transpose[m]+m
]
]*)


(*Position[Map[Tag,base2,{2}], Tag@{0,0,1}][[1,1]]*)
