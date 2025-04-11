#!/usr/bin/env wolframscript
(* ::Package:: *)

Import["../../Mathematica_packages/QMB.wl"]


PositiveParitySubspaceBasis[L_]:=DeleteDuplicatesBy[Tuples[{0,1},L],First[Sort[{#,Reverse[#]}]]&]


time1=AbsoluteTime[];


L = 16;

base = PositiveParitySubspaceBasis[L];
tagsbase = Tag /@ base;

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
];


time2=AbsoluteTime[];


Print["Construcci\[OAcute]n de Hx tom\[OAcute] "<>ToString[NumberForm[time2-time1,{Infinity,2}]]<>" s"];


Hz=DiagonalMatrix[Total[(-1)^#]&/@base, TargetStructure->"Sparse"];
JJ=DiagonalMatrix[Total[(-1)^Differences[#]]&/@base, TargetStructure->"Sparse"];
Hamiltonian[hx_,hz_,J_]:=hx*Hx+hz*Hz-J*JJ;


points = Tuples[Join[Flatten[Table[Subdivide[#[[i]],#[[i+1]],3][[1;;3]],{i,9}]],{2.95}],2]&[Subdivide[0.05,2.95,9]];

time1=AbsoluteTime[];

LaunchKernels[56];
Print["Launched kernels : ", Length[ParallelKernels[]]];
eigenvalues = 
	ParallelTable[
		Module[{hz, J, H},
		{hz, J} = i;
		H = Normal[Hamiltonian[1., hz, J]];
		Flatten[{{hz, J}, Eigenvalues[H]}]
		]
	, {i, points},
	DistributedContexts->Full, Method->"CoarsestGrained"];
CloseKernels[];

time2=AbsoluteTime[];
Print["Tiempo diagonalizaciones: "<>ToString[NumberForm[time2-time1,{Infinity,2}]]<>" s"];

Export["data_eigenvalues_wisniacki_L_"<>ToString[L]<>"_hx_1.csv",Prepend[eigenvalues,{"hz,J,E_1,E_2,...,E_{2^L}"}],"CSV"];


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
