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
(*Load right-regular representation of generators*)


Adata = Import[fprefix<>"_A.reg","Table"];
Bdata = Import[fprefix<>"_B.reg","Table"];
iAdata = Import[fprefix<>"_iA.reg","Table"];
iBdata = Import[fprefix<>"_iB.reg","Table"];

\[CapitalPi]["A"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Adata, Null], {dim,dim} ] ;
\[CapitalPi]["B"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Bdata, Null], {dim,dim} ] ;
\[CapitalPi]["iA"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iAdata, Null], {dim,dim} ] ;
\[CapitalPi]["iB"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iBdata, Null], {dim,dim} ] ;


Agraph = DirectedEdge[#[[1]],#[[2]]] &/@  \[CapitalPi]["A"]["ExplicitPositions"];
Bgraph = DirectedEdge[#[[1]],#[[2]]] &/@  \[CapitalPi]["B"]["ExplicitPositions"];
iAgraph = DirectedEdge[#[[1]],#[[2]]] &/@  \[CapitalPi]["iA"]["ExplicitPositions"];
iBgraph = DirectedEdge[#[[1]],#[[2]]] &/@  \[CapitalPi]["iB"]["ExplicitPositions"];


gr = Graph3D[Agraph~Join~Bgraph, GraphLayout->{ "SpringElectricalEmbedding" ,"RepulsiveForcePower"-> -2, "SpringConstant"-> 10}]


grphpl[\[Phi]_]:=GraphPlot3D[gr,ViewPoint->3.4{Cos[\[Phi]],Sin[\[Phi]],0},SphericalRegion->True , ImageSize -> {500, 500}]


(*phis= Subdivide[0,\[Pi]/4, 40];
phis = phis ~ Join ~ Drop[ Reverse[phis], 1 ]*)


(*plts = grphpl[#]&/@ phis*)


(* Export["anim.gif",plts, AnimationRate-> 25] *)


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


PDInitCanvas[];
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


Ageo = PDGeodesic[zs[[#[[1]]]],zs[[#[[2]]]],20] &/@  \[CapitalPi]["A"]["ExplicitPositions"];
Bgeo = PDGeodesic[zs[[#[[1]]]],zs[[#[[2]]]],20] &/@  \[CapitalPi]["B"]["ExplicitPositions"];


PDInitCanvas[];
PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo; 
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


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
