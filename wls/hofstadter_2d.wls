#!/usr/bin/env wolframscript
(* ::Package:: *)

n = 2;
p = 2;
mod = p^n;
index = Range[0,mod-1];
hilbertspace= Tuples[index,2];
d = Length[hilbertspace];


\[Gamma][1] = SparseArray[ { # , FirstPosition[ hilbertspace, Mod[ hilbertspace[[#]]+ {1,0} , mod ]][[1]]} ->1  &/@ Range[d] , {d,d}]


\[Gamma][2] = SparseArray[ { # , FirstPosition[ hilbertspace, Mod[ hilbertspace[[#]]+ {0,1} , mod ]][[1]]} ->1  &/@ Range[d] , {d,d}]


H = -(\[Gamma][1] + \[Gamma][2]) ;
H =H + ConjugateTranspose[H] //N //Normal ;


spec = Eigenvalues[H];


MagSpec[n_, k_]:=(

	s = Exp[n (2 \[Pi] I)/k];
	
	\[CapitalSigma][1] := If[k>1,SparseArray[{Band[{2,1}]->1},{k,k}], SparseArray[{},{k,k}]]  + SparseArray[{ {1,k}-> 1}, {k,k}];
	\[CapitalSigma][2] := SparseArray[ {#,#}->s^(#-1) &/@ Range[k], {k,k}];
	
	H = -Sum[ KroneckerProduct[\[CapitalSigma][i], \[Gamma][i]],{i,1,2}] ;
	
	H = H + ConjugateTranspose[H] //N //Normal ;
	
	spec = {n/k,#} &/@ Eigenvalues[H];
	
	Return[spec];
)


k = 101;
Print[d k];


PrimeQ[83]


specs = MagSpec[#, k]&/@Range[k] ;


ListPlot[ArrayFlatten[specs,1],AspectRatio->1,Frame->True,PlotLabel->"g=1, \!\(\*SuperscriptBox[\(p\), \(n\)]\)=4, k=101", FrameLabel->{"Flux","Energy"}]
