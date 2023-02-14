#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Generating the ring extension*)


(* ::Subtitle:: *)
(*Iterative algorithm*)


(* ::Text:: *)
(*./wls/generate_ring_reduction.wls*)
(**)
(*Author: 	Fabian R. Lux*)
(*Date:	2/14/2023*)
(**)
(*This Wolfram language script can be used to construct the irreducible minimal polynomials Subscript[\[CapitalPsi], n][x] which have Subscript[\[Xi], n ]= 2 cos(2\[Pi]/n) as a root. *)
(*Subscript[\[CapitalPsi], n][x]  has degree d=\[Phi](n)/2, where \[Phi](n) is Euler's totient function. Further, Subscript[\[CapitalPsi], n][x]  is monic and has integer coefficients.*)
(*This means we can express Subscript[\[Xi], n]^d as a polynomial of degree <d with integer coefficients. *)
(**)
(*The output of this script is a file "./wls/ring_reduction.inp" which contains a list of these polynomials for 2<=n<=250.*)
(**)
(*The script can also be used to generate a L AT EX table of the irreducible polynomials (last section of this script).*)
(**)
(**)


nmax := 250;


P[n_]:=2 ChebyshevT[n, x/2]


Clear[\[CapitalPsi]];
\[CapitalPsi][n_]:=\[CapitalPsi][n]=(
	div = Drop[Divisors[n],-1];
	
	rhs = If[EvenQ[n], 
	P[n/2+1]-P[n/2-1],  
	P[(n-1)/2+1]-P[(n-1)/2]];
	
	recursion = Times@@(\[CapitalPsi][#]&/@ div );
	
	Return[Expand[ FullSimplify[ rhs/recursion ]]]
);


\[CapitalPsi][#]&/@ Range[nmax];


(* ::Text:: *)
(*The following function gives the representation of Subscript[\[Xi], n]^d in terms of the promised polynomial with degree less than d.*)


r[n_]:= {n, EulerPhi[n]/2} ~ Join~PadRight[CoefficientList[Expand[ x^(EulerPhi[n]/2)-\[CapitalPsi][n] ],x],EulerPhi[n]/2]


prefix = "IRREDUCIBLE POLYNOMIALS\n\n";
prefix = prefix <> "\tAuthor: Fabian R. Lux\n";
prefix = prefix <>"\tDate: 2/13/2023\n\n";
prefix = prefix <> "\tLet Psi(n,x) represent the irreducible, minimal polynomial in Z[x] which has xi(n)=2 cos(2 pi / n) as a root.\n";
prefix = prefix <> "\tThen Psi(n,x) has degree d=phi(n)/2, where phi(n) is Euler's totient function.\n\n";
prefix = prefix <> "\tThis means we can express xi(n)^d as a polynomial xi(n)^d = a(0)*xi(n)^0 + a(1)*xi(n)^1  + ... + a(d-1)*xi(n)^(d-1).\n\n";
prefix = prefix <> "\tThe following list is formated such that each line starts off by giving the index n followed by the degree d. \n";
prefix = prefix <> "\tThen follows a sequence of integers which correspond to the coefficients a(0), a(1), ..., a(d-1).\n\n";
prefix = prefix <> "BEGIN\n";


postfix ="\nEND";


rows = r[#] &/@ Range[3,nmax];


FormatRow[row_]:=StringRiffle[ ToString[#]&/@row," "] 


lines = StringRiffle[ FormatRow[# ]&/@ rows, "\n" ] ;


data = prefix <> lines <>postfix ;


dir = If[$Notebooks, NotebookDirectory[], Directory[] ];
SetDirectory[dir];
Export["./ring_reduction.inp", data, "Text"]


(* ::Chapter:: *)
(*T EX Table *)


pqlist = DeleteCases[ ArrayFlatten[ Table[ { {p,q}, 1/p+1/q+1/2<1}, {p,1,20},{q,1,20}],1], {{p_,q_}, False}] ;


pqcurated= SortBy[DeleteCases[ ArrayFlatten[ Table[ { p,q,  \[Pi]-\[Pi]/p-\[Pi]/q-\[Pi]/2, \[CapitalPhi][2 p q]}, {p,1,8},{q,1,p}], 1], {p_,q_,x_/;(x<=0),y_}],#[[3]]&] ;
Length[pqcurated]


pqcurated


TeXPrefix ="\\begin{table*}[t]
\\caption{ {\\bf Irreducible polynomials for $\\xi_{2pq}=2 \\cos(\\pi/(pq))$}. For $p\\leq 8$ and $q \\leq p$, we list all admissible hyperbolic triangles, sorted by ascending hyperbolic area.
The last column contains the minimal irreducible polynomial $\\Psi_{2pq}(x) $, which has $\\xi_{2pq}$ as a root.}
\\centering
\\scalebox{0.8}{
\\begin{tabular}{cccl}
\\toprule
$p$  & $q$ & Area & $\\Psi_{2pq}(x) $  \\\\
\\midrule
";


TeXPostFix="
\\bottomrule
\\end{tabular}}
\\label{table:irreducible_polynomials}
\\end{table*}";


FormatRow[row_]:=StringRiffle[ "$"<>ToString[ TeXForm[#]]<>"$" &/@row, " & " ] <> " \\\\ \n"


TeXTable = TeXPrefix <> StringJoin[ FormatRow[#]&/@pqcurated ] <> TeXPostFix



