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


dir = If[$Notebooks, NotebookDirectory[], Directory[] ];
SetDirectory[dir];


pmax= 20;
qmax = 20;


P[n_]:=2 ChebyshevT[n, x/2]


pqlist = #[[1]]&/@DeleteCases[ ArrayFlatten[ Table[ { {p,q}, 1/p+1/q+1/2<1}, {p,1,pmax},{q,1,qmax}],1], {{p_,q_}, False}] ;


Clear[n]
n[p_,q_]:=(
	res = 2 p q  ;
	
	If[ p==q , res = 2p; ];
	
	If[ q==3 , res = 2p ;  ];
	
	If[ p==3 , res =2 q ; ];
	
	Return[res];
)


ExpandPoly[poly_,p_, q_]:=PadRight[CoefficientList[Expand[ poly ],x],EulerPhi[n[p,q]]/2]


X[p_,q_]:=(
	
	res = {{-1, P[q],0},{0, 1, 0 },{0,0,1}}; 
	
	If[ p==q , res = {{-1, x,0},{0, 1, 0 },{0,0,1}};  ];
	
	If[ q==3 , res = {{-1, x,0},{0, 1, 0 },{0,0,1}};  ];
	
	If[ p==3 , res = {{-1, 1,0},{0, 1, 0 },{0,0,1}};  ];

	Return[res];
)

Y[p_,q_]:=(
	
	res = {{1, 0,0},{P[q], -1, P[p] },{0,0,1}};
	
	If[ p==q , res = {{1, 0,0},{x, -1, x },{0,0,1}}; ];
	
	If[ q==3 , res = {{1, 0,0},{x, -1, 1 },{0,0,1}}; ];
	
	If[ p==3 , res = {{1, 0,0},{1, -1, x },{0,0,1}};];

	Return[res];
)

Z[p_,q_]:=(
	
	res = {{1, 0,0},{0, 1, 0 },{0,P[p],-1}};
	
	If[ p==q , res = {{1, 0,0},{0, 1, 0 },{0,x,-1}};];
	
	If[ q==3 , res = {{1, 0,0},{0, 1, 0 },{0,1,-1}}; ];
	
	If[ p==3 , res = {{1, 0,0},{0, 1, 0 },{0,x,-1}};];

	Return[res];
)


XFormat[p_,q_]:= ExpandPoly[#, p, q] &/@ Flatten[ Expand[ X[p,q] ]]
YFormat[p_,q_]:= ExpandPoly[#, p, q] &/@ Flatten[ Expand[ Y[p,q] ]]
ZFormat[p_,q_]:= ExpandPoly[#, p, q] &/@ Flatten[ Expand[ Z[p,q] ]]


XExport[{p_,q_}]:=Export["./generators/"<>ToString[p]<>"_"<>ToString[q]<>".X", XFormat[p,q], "Table", "FieldSeparators" -> " "]
YExport[{p_,q_}]:=Export["./generators/"<>ToString[p]<>"_"<>ToString[q]<>".Y", YFormat[p,q], "Table", "FieldSeparators" -> " "]
ZExport[{p_,q_}]:=Export["./generators/"<>ToString[p]<>"_"<>ToString[q]<>".Z", ZFormat[p,q], "Table", "FieldSeparators" -> " "]


XExport[#]&/@ pqlist;
YExport[#]&/@ pqlist;
ZExport[#]&/@ pqlist;
