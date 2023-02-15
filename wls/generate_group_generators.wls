#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Generating the group generators*)


(* ::Subtitle:: *)
(*Using Chebyshev polynomials*)


(* ::Text:: *)
(*./wls/generate_group_generators.wls*)
(**)
(*Author: 	Fabian R. Lux*)
(*Date:	2/14/2023*)
(**)
(*Defines and export the group generators for all admissible p and q for p,q<20. The result is stored in the folder ./wls/generators*)
(**)
(*The X generators are found in <p>_<q>.X*)
(*The Y generators are found in <p>_<q>.Y*)
(*The Z generators are found in <p>_<q>.Z*)
(**)
(**)


pmax= 20;
qmax = 20;


P[n_]:=2 ChebyshevT[n, x/2]


pqlist = #[[1]]&/@DeleteCases[ ArrayFlatten[ Table[ { {p,q}, 1/p+1/q+1/2<1}, {p,1,pmax},{q,1,qmax}],1], {{p_,q_}, False}] ;


ExpandPoly[poly_,n_]:=PadRight[CoefficientList[Expand[ poly ],x],EulerPhi[n]/2]


X[p_,q_]:={{-1, P[q],0},{0, 1, 0 },{0,0,1}};
Y[p_,q_]:={{1, 0,0},{P[q], -1, P[p] },{0,0,1}};
Z[p_,q_]:={{1, 0,0},{0, 1, 0 },{0,P[p],-1}};


XFormat[p_,q_]:= ExpandPoly[#,2 p q] &/@ Flatten[ Expand[ X[p,q] ]]
YFormat[p_,q_]:= ExpandPoly[#,2 p q] &/@ Flatten[ Expand[ Y[p,q] ]]
ZFormat[p_,q_]:= ExpandPoly[#,2 p q] &/@ Flatten[ Expand[ Z[p,q] ]]


dir = If[$Notebooks, NotebookDirectory[], Directory[] ];
SetDirectory[dir];


XExport[{p_,q_}]:=Export["./generators/"<>ToString[p]<>"_"<>ToString[q]<>".X", XFormat[p,q], "Table", "FieldSeparators" -> " "]
YExport[{p_,q_}]:=Export["./generators/"<>ToString[p]<>"_"<>ToString[q]<>".Y", YFormat[p,q], "Table", "FieldSeparators" -> " "]
ZExport[{p_,q_}]:=Export["./generators/"<>ToString[p]<>"_"<>ToString[q]<>".Z", ZFormat[p,q], "Table", "FieldSeparators" -> " "]


XExport[#]&/@ pqlist;
YExport[#]&/@ pqlist;
ZExport[#]&/@ pqlist;
