#!/usr/bin/env wolframscript
(* ::Package:: *)

Get["~/Documents/libs/QMB.wl"];
Get["/home/jadeleon/Documents/chaos_meets_channels/Mathematica_packages/Chaometer.wl"]
LogSpace[start_,end_,n_]:=Exp@Subdivide[Log[start],Log[end],n-1];


(*points = Tuples[{#,#}]&[Subdivide[0.05, 2.95, 9]];*)
points = Tuples[Subdivide[0.05,5,29], 2];
lenpoints = Length[points];


t0=0.1;tf=1000;n=100(* cantidad de tiempos *);
t=LogSpace[t0,tf,100];


L=7;
Jz=1.;\[Omega]=0.;d=3;


Do[
{\[Epsilon]d,Jxy} = points[[i]];
(*If[FileExistsQ["~/Documents/chaos_meets_channels/ja_files/data/chois/xxz/L_7_Jz_1_omega_0_epsilon_"<>ToString[NumberForm[\[Epsilon]d,{Infinity,2}]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>"/100.mat"], Continue[]];*)

t1 = AbsoluteTime[];
	Do[
	H=LeaSpinChainHamiltonian[Jxy,Jz,\[Omega],\[Epsilon]d,L,d];
	{eigenvals,eigenvecs}=Chop[Eigensystem[H]];
	
	x=RandomInteger[100000];SeedRandom[x];
	\[Psi]=RandomChainProductState[L-1];
	A=KroneckerProduct[IdentityMatrix[2^(L-1),SparseArray]/2,Dyad[\[Psi]]];
	
	P=Transpose[eigenvecs];Pinv=Conjugate[eigenvecs];
	ClearAll[U];U[t_?NumericQ]:=P . DiagonalMatrix[Exp[-I*eigenvals*t],TargetStructure->"Sparse"] . Pinv;
	
	chois=Table[ChoiMatrix[U[j],A,L],{j,t}];
	Export["~/Documents/chaos_meets_channels/ja_files/data/chois/xxz/L_7_Jz_1_omega_0_epsilon_"<>ToString[NumberForm[\[Epsilon]d,{Infinity,2}]]<>"_Jxy_"<>ToString[NumberForm[Jxy,{Infinity,2}]]<>"/"<>StringPadLeft[ToString[k],3,"0"]<>".mat",
	Prepend[Flatten/@chois,{"Choi matrices from t=0 to t=1000 log spaced, 100 times. Initial state of environment generated with SeedRandom["<>ToString[x]<>"];\[Psi]=RandomChainProductState[L-1];"}],"MAT"];
	, {k,100}];
t2 = AbsoluteTime[]; Print["(hz, J) = (", ToString[NumberForm[hz, {Infinity, 2}]], ", ", ToString[NumberForm[J, {Infinity, 2}]], ") i = ", i, ": ", ToString[NumberForm[t2-t1, {Infinity, 1}]], "s"];

, {i, lenpoints}]
