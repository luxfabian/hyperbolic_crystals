#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*This script is only for the {8,3} triangle group*)


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
(*Load right-regular representation of generators*)


Adata = Import[fprefix<>"_A.reg","Table"];
Bdata = Import[fprefix<>"_B.reg","Table"];
iAdata = Import[fprefix<>"_iA.reg","Table"];
iBdata = Import[fprefix<>"_iB.reg","Table"];

\[CapitalPi]["A"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Adata, Null], {dim,dim} ] ;
\[CapitalPi]["B"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Bdata, Null], {dim,dim} ] ;
\[CapitalPi]["iA"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iAdata, Null], {dim,dim} ] ;
\[CapitalPi]["iB"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iBdata, Null], {dim,dim} ] ;


ParseWord[word_]:=(
z = z0;
( z = Moebius[g[#],z])&/@ Reverse[ Characters[ word ] ];
Return[z]
)


WordForOrbit = 3;

\[Gamma]1word = "BAAAA";
\[Gamma]1segments = \[CapitalPi][#]&/@   Reverse[ Characters[ \[Gamma]1word ] ] ;

Clear[\[Gamma]path];
\[Gamma]1len = Length[Reverse[ Characters[ \[Gamma]1word ] ]];
\[Gamma]path[0]:= SparseArray[ IdentityMatrix[dim] ];
\[Gamma]path[i_]:= \[Gamma]path[i]= \[Gamma]1segments[[i]] . \[Gamma]path[i-1];


Print["Chosen word: ", lexikon[[WordForOrbit]]];

wordvec = Table[ KroneckerDelta[WordForOrbit, i],{i,1,dim}];


orbit = Table[ (\[Gamma]path[i] . wordvec) . zs, {i,0,\[Gamma]1len}];

geo=PDGeodesic[orbit[[#]],orbit[[ #+1]],50, {Thick,fuchscolor}]&/@ Range[Length[orbit]-1];

PDInitCanvas[];
orbitpts = PDHyperbolicPoint[#]&/@orbit;

PDInitCanvas[];
PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo; 
PDUpdateCanvas[#]&/@pts;
PDUpdateCanvas[#]&/@ fuchsiangeo;
gfx=PDShow[]


Length[zs]





Length][


\[Gamma]10





(* ::Input:: *)
(*\[CapitalPi][#]&/@{"A","A","A","A","B","A","A","A","A\","A ","B ","A "}*)


DirectedEdge[#[[2]],#[[1]]] &/@  \[CapitalPi]["A"]["ExplicitPositions"];


(* ::Text:: *)
(*Setup the generators*)


\[Alpha] =2 \[Pi]/p;
\[Beta]=2 \[Pi]/q;
r0  = Sqrt[Cos[(\[Alpha]+\[Beta])/2]/Cos[(\[Alpha]-\[Beta])/2]];

g["A"] =  { { Exp[I \[Alpha] /2],0},{0, Exp[-I \[Alpha]/2]}} ;
g["B"] = 1/(1-r0^2) { { Exp[I \[Beta] /2]-r0^2 Exp[-I \[Beta]/2],r0 (1-Exp[I \[Beta]]) Exp[I (\[Alpha]-\[Beta])/2]},{r0 (1-Exp[-I \[Beta]]) Exp[-I (\[Alpha]-\[Beta])/2],Exp[-I \[Beta]/2]-r0^2 Exp[I \[Beta]/2]}} ;


\[Alpha]f =2 \[Pi]/8;
\[Beta]f=2 \[Pi]/8;
r0f  = Sqrt[Cos[(\[Alpha]f+\[Beta]f)/2]/Cos[(\[Alpha]f-\[Beta]f)/2]];

gf["A"] =  { { Exp[I \[Alpha]f /2],0},{0, Exp[-I \[Alpha]f/2]}} ;
gf["B"] = 1/(1-r0f^2) { { Exp[I \[Beta]f /2]-r0f^2 Exp[-I \[Beta]f/2],r0 (1-Exp[I \[Beta]f]) Exp[I (\[Alpha]f-\[Beta]f)/2]},{r0f (1-Exp[-I \[Beta]f]) Exp[-I (\[Alpha]f-\[Beta]f)/2],Exp[-I \[Beta]f/2]-r0f^2 Exp[I \[Beta]f/2]}} ;


(* ::Text:: *)
(*Group elements act by Moebius transformations*)


Moebius[ {{a_,b_},{c_,d_}}, z_]:= (a z + b )/(c z +d )


ParseWord[word_]:=(
z = z0;
( z = Moebius[g[#],z])&/@ Reverse[ Characters[ word ] ];
Return[z]
)


z0 =  Exp[I (\[Alpha]/2 + 0 \[Pi]/6)](r0 + 0.07 Exp[I \[Pi]]);

zs= N[ParseWord[#]]&/@lexikon;
pts = PDHyperbolicPoint[#]&/@zs;

PDInitCanvas[];
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


\[Gamma]1=MatrixPower[ g["A"] . g["B"] . MatrixPower[g["A"],4], 2]//Simplify;
\[Gamma][n_]:=MatrixPower[ g["A"], n-1] . \[Gamma]1 . MatrixPower[ Inverse[g["A"]], n-1]





loop = \[Gamma][#] &/@ Range[8] //N;


loop = {\[Gamma][1], Inverse[\[Gamma][2]], \[Gamma][3], Inverse[\[Gamma][4]], Inverse[\[Gamma][1]],\[Gamma][2], Inverse[\[Gamma][3]],\[Gamma][4] } //N ;


Dot@@loop// Chop


Clear[\[Gamma]path];


seed = 0;
\[Gamma]path[0]:= IdentityMatrix[2];
\[Gamma]path[i_]:= \[Gamma]path[i]= \[Gamma]path[i-1] . Inverse[ Reverse[loop][[i]] ];



orbit = Table[ Moebius[\[Gamma]path[i],seed], {i,0,Length[loop]}];

geo=PDGeodesic[orbit[[#]],orbit[[ #+1]],300, {Thick,Black}]&/@ Range[Length[orbit]-1];

PDInitCanvas[];
orbitpts = PDHyperbolicPoint[#]&/@orbit;

PDInitCanvas[];
PDUpdateCanvas[#]&/@geo;
PDUpdateCanvas[#]&/@ fuchsiangeo;
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


fuchsian= Exp[I \[Alpha]f #] Exp[I \[Alpha]f/2]r0f &/@ Range[8]//N ;


fuchscolor=RGBColor["#9f0a0a"];


fuchsiangeo=PDGeodesic[fuchsian[[#]],fuchsian[[ Mod[#,8]+1]],50, {Thick,fuchscolor}]&/@ Range[8];


PDInitCanvas[];
PDUpdateCanvas[#]&/@ fuchsiangeo;
gfx=PDShow[]


PDInitCanvas[];
pts = PDHyperbolicPoint[#]&/@zs;
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


Ageo = PDGeodesic[zs[[#[[1]]]],zs[[#[[2]]]],20] &/@  \[CapitalPi]["A"]["ExplicitPositions"];
Bgeo = PDGeodesic[zs[[#[[1]]]],zs[[#[[2]]]],20] &/@  \[CapitalPi]["B"]["ExplicitPositions"];


PDInitCanvas[];
PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo; 
PDUpdateCanvas[#]&/@pts;
PDUpdateCanvas[#]&/@ fuchsiangeo;
gfx=PDShow[]


Export["lattice.png",gfx, ImageResolution -> 300, Background -> None]


PartialPlot[n_]:=(
PDInitCanvas[];
PDSetVertexColor[Black];
PDUpdateCanvas[#]&/@( pts[[1;;n]]);
Return[PDShow[]];
);


(* ::Input:: *)
(**)


allplots = PartialPlot[#]&/@ Range[Length[lexikon]] ;


Export["test.gif",allplots]


vI


col[z_]:= makeBlend[""][Rescale[Arg[z],{-\[Pi],\[Pi]}]] 


Get["https://pastebin.com/raw/qcyE3vzF"]


[0]


PDVertexColor


sc = 1;


bl[x_]:=Blend[{White,RGBColor["#990000"]}, Rescale[x,{Log10[1],Log10[1+sc]}]]


bl[x_]:=Blend[{White,RGBColor["#990000"]}, x]


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


ldos = Flatten[ Import["../py/ldos_abs_mag.dat", "Data"] ];

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


Export["ldos_crystal_0.2.png", gfx2, ImageResolution-> 400]


blgfx=BarLegend[{
bl,
  {0,1}
  },LegendLayout->"Column"]


Export["ldos_color_bar.png", blgfx, ImageResolution-> 300]
