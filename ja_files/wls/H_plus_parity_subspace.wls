#!/usr/bin/env wolframscript
(* ::Package:: *)

PositiveParitySubspaceBasis[L_]:=DeleteDuplicatesBy[Tuples[{0,1},L],First[Sort[{#,Reverse[#]}]]&]


time1=AbsoluteTime[];


L=17;
base=PositiveParitySubspaceBasis[L];

(*Esto seguro se puede optimizar, pero queda para despu\[EAcute]s*)
Hx=Module[{m,dimSubspace=Length[base],q,l},
m=SparseArray[Prepend[
Table[
tmp1=Select[Join[ConstantArray[j,L]-IdentityMatrix[L],If[PalindromeQ[j],{},ConstantArray[Reverse[j],L]-IdentityMatrix[L]]],And@@NonNegative[#]&];
q=ReplaceAll[PalindromeQ[j],{True->0,False->1}];
Total[
Table[
l=ReplaceAll[PalindromeQ[k],{True->0,False->1}];
SparseArray[Flatten[Position[base,If[MemberQ[base,k],k,Reverse[k]]]]->2^(-(q+l)/2),dimSubspace]
,{k,tmp1}]
]
,{j,base[[2;;]]}]
,SparseArray[{},dimSubspace]
]];
ConjugateTranspose[m]+m
];


time2=AbsoluteTime[];


Print["Construcci\[OAcute]n de Hx tom\[OAcute] ", NumberForm[time2-time1,{Infinity,2}], " s"];


Hz=DiagonalMatrix[Total[(-1)^#]&/@base,TargetStructure->"Sparse"];


JJ=DiagonalMatrix[Total[(-1)^Differences[#]]&/@base,TargetStructure->"Sparse"];


Hamiltonian[hx_,hz_,J_]:=hx*Hx+hz*Hz-J*JJ;


sizeH=ByteCount[H = Hamiltonian[1., 0.5, 1.]];


Print["Tama\[NTilde]o en memoria de H: ", sizeH/10.^9, " GB"];


time1=AbsoluteTime[];


memoryDiagonalization = MaxMemoryUsed[Eigenvalues[H]];


time2=AbsoluteTime[];


Print["Memoria para diagonalizar: ", memoryDiagonalization/10.^9, " GB"];


Print["Tiempo para diagonalizar: ", NumberForm[time2-time1,{Infinity,2}], " s"];
