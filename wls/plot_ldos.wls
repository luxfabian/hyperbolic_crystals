#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[NotebookDirectory[]];

Import["poincare.wls"];


(* ::Text:: *)
(*Get info about group*)


{p,q,n} = #[[2]]&/@Import["../data/group_specs.inp", "Table"];

If[n>0,
fprefix = "../data/"<>ToString[p]<>"_"<>ToString[q]<>"_modulo_"<>ToString[n],
fprefix = "../data/"<>ToString[p]<>"_"<>ToString[q]<>"_open_"<>ToString[-n]];

{dim}=Import[fprefix<>".info", "Data"]//ToExpression;

Print["Dimension of Hilbert space: ", dim]


(* ::Text:: *)
(*Read the basis file*)


lexikon = Import[fprefix<>".words", "Data"];
lexikon =#[[2]]&/@ SortBy[lexikon,First];


(* ::Text:: *)
(*Setup the generators*)


\[Alpha] =2 \[Pi]/p;
\[Beta]=2 \[Pi]/q;
r0  = Sqrt[Cos[(\[Alpha]+\[Beta])/2]/Cos[(\[Alpha]-\[Beta])/2]];

g["A"] =  { { Exp[I \[Alpha] /2],0},{0, Exp[-I \[Alpha]/2]}} ;
g["B"] = 1/(1-r0^2) { { Exp[I \[Beta] /2]-r0^2 Exp[-I \[Beta]/2],r0 (1-Exp[I \[Beta]]) Exp[I (\[Alpha]-\[Beta])/2]},{r0 (1-Exp[-I \[Beta]]) Exp[-I (\[Alpha]-\[Beta])/2],Exp[-I \[Beta]/2]-r0^2 Exp[I \[Beta]/2]}} ;


(* ::Text:: *)
(*Seed coordinate*)


z0 =  Exp[I \[Alpha]/2](r0 + 0.12 Exp[I \[Pi]]);


(* ::Text:: *)
(*Group elements act by Moebius transformations*)


Moebius[ {{a_,b_},{c_,d_}}, z_]:= (a z + b )/(c z +d )


ParseWord[word_]:=(
z = z0;
( z = Moebius[g[#],z])&/@ Reverse[ Characters[ word ] ];
Return[z]
)


zs= N[ParseWord[#]]&/@lexikon;
pts = PDHyperbolicPoint[#]&/@zs;


PDHyperbolicPoint2[zi_,vi_,vmax_]:=(
	Clear[d];
	d = ArcTanh[Abs[zi]];
	
	color = RGBColor["#990000"];
	
	ptsize=0.07;
	
	radius = ptsize (*  Exp[-d]; *);
	
	If[ Abs[vi/vmax] > -0.001,
	point = Graphics[{color, Opacity[ Abs[vi/vmax]], Disk[{Re[zi],Im[zi]},radius]}];
	
	Return[point];
	];
);


ldos = Flatten[ Import["../py/out/ldos_abs_mag.dat", "Data"] ];

vmax = Max[ Abs[ldos] ];

pts = DeleteCases[PDHyperbolicPoint2[ zs[[#]], ldos[[#]], vmax] &/@ Range[dim], Null];


PDInnerColorOfPoincareDisk:=White;
PDOuterColorOfPoincareDisk:=White;

PDInitCanvas[];
PDVertexColor=RGBColor[  "#1e81b0" ];
PDUpdateCanvas[PDHyperbolicPoint[#]]&/@zs;
PDUpdateCanvas[#]&/@pts;
PDUpdateCanvas[ PDGeodesic[ 0, Exp[I \[Alpha]/2],10] ]
PDUpdateCanvas[ PDGeodesic[ 0, Exp[I \[Alpha]/2 + I 2\[Pi]/3],10] ]
PDUpdateCanvas[ PDGeodesic[ 0, Exp[I \[Alpha]/2+ I 2 2\[Pi]/3],10] ]

gfx=PDShow[];
gfx2=Show[gfx,ImageSize->100]


Export["./out/ldos_crystal_0.2.png", gfx2, ImageResolution-> 400]


blgfx=BarLegend[{
bl,
  {0,1}
  },LegendLayout->"Column"]


bl[x_]:=Blend[{White,RGBColor["#990000"]}, x]
Export["./out/ldos_color_bar.png", blgfx, ImageResolution-> 300]
