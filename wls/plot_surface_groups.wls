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


z0 =  Exp[I \[Alpha]/2]r0;


(* ::Text:: *)
(**)
(*Group elements act by Moebius transformations*)


Moebius[ {{a_,b_},{c_,d_}}, z_]:= (a z + b )/(c z +d )


ParseWord[word_]:=(
z = z0;
( z = Moebius[N[g[#]],z])&/@ Reverse[ Characters[ word ] ];
Return[z]
)


ParseDualWord[word_]:=(
z = 0;
( z = Moebius[N[g[#]],z])&/@ Reverse[ Characters[ word ] ];
Return[z]
)


zs= N[ParseWord[#]]&/@lexikon;
pts = PDHyperbolicPoint[#]&/@zs;

zsdual= N[ParseDualWord[#]]&/@lexikon;
ptsdual = PDHyperbolicPoint[#]&/@zsdual;


PDInitCanvas[];
PDUpdateCanvas[#]&/@ptsdual;
gfx=PDShow[]


Ageo = PDGeodesic[zs[[#[[1]]]],zs[[#[[2]]]],20] &/@  \[CapitalPi]["A"]["ExplicitPositions"];
Bgeo = PDGeodesic[zs[[#[[1]]]],zs[[#[[2]]]],20] &/@  \[CapitalPi]["B"]["ExplicitPositions"];

Ageodual = PDGeodesic[zsdual[[#[[1]]]],zsdual[[#[[2]]]],20,{Thin,PDVertexColor}] &/@  \[CapitalPi]["A"]["ExplicitPositions"];
Bgeodual = PDGeodesic[zsdual[[#[[1]]]],zsdual[[#[[2]]]],20,{Thin,PDVertexColor}] &/@  \[CapitalPi]["B"]["ExplicitPositions"];


PDInitCanvas[];
PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo;
(*PDUpdateCanvas[#]&/@pts;*)

PDUpdateCanvas[#]&/@Ageodual;
PDUpdateCanvas[#]&/@Bgeodual;


fd = {0, z0, GeodesicMidpoint[z0, Moebius[g["A"],z0]] };
ptsfd = PDHyperbolicPoint[#]&/@fd;

style=Black;
geofd = { PDGeo[fd[[1]], fd[[2]],30,style],PDGeo[fd[[2]], fd[[3]],30,style],PDGeo[fd[[3]], fd[[1]],30,style] };

PDUpdateCanvas[#]&/@geofd;

PDUpdateCanvas[#]&/@ptsdual;

gfx=PDShow[]


Export["octogon.png",gfx,ImageResolution-> 300]


z = z0;


polygon = {z0};
Iterate[]:= (
	polygon = polygon ~ Join ~ { Moebius[g["A"],Last[polygon]]};
);
Iterate[]&/@ Range[p-1];
polygonpts = PDHyperbolicPoint[#]&/@polygon;
polygongeos = Table[ PDGeodesic[polygon[[i]],polygon[[i+1]],20] , {i,1,7}]~ Join ~{PDGeodesic[polygon[[8]],polygon[[1]],20]};












Length[polygon]


polygon





Mod[2,2]





PDInitCanvas[];
PDUpdateCanvas[#]&/@geos;
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


polygonT = Moebius[Inverse[g["A"]] . Inverse[g["A"]] . g["A"] . Inverse[ g["B"] ] . g["A"] . g["A"],#] &/@ polygon //N;
geosT = Table[ PDGeodesic[polygonT[[i]],polygonT[[i+1]],20] , {i,1,7}]~ Join ~{PDGeodesic[polygonT[[8]],polygonT[[1]],20]};
polygonT





\[Gamma]1 = g["B"] . g["A"] . MatrixPower[g["A"],p/2 ] ;


\[Gamma]1 = g["A"] . g["B"] . MatrixPower[g["A"],p/2 ] ;


\[Gamma]1B = g["B"] . g["A"] . MatrixPower[g["B"],p/2 ] ;


\[Gamma][n_]:= MatrixPower[g["A"], n-1 ] . \[Gamma]1 . MatrixPower[g["A"], -n+1 ];


\[Gamma]B[n_]:= MatrixPower[g["B"], n-1 ] . \[Gamma]1B . MatrixPower[g["B"], -n+1 ];


polygonT


Solve[ Moebius[g["A"],zx + I zy]==0, {zx,zy}]


Ahlfors[z1_,z2_] := Sqrt[(1 - Abs[z1]^2)(1 - Abs[z2]^2) + Abs[z1-z2]^2]


GeodesicMidpoint[z1_, z2_]:=(z2 (1- Abs[z1]^2)+z1 (1- Abs[z2]^2))/(1- Abs[z1]^2  Abs[z2]^2+Ahlfors[z1,z2] Sqrt[(1-Abs[z1]^2)(1-Abs[z2]^2)])


(* ::Text:: *)
(*Modified geodesic routine*)


PDGeo[x_, y_, pts_,style_]:=(

	(*Translate point x to 0 (then the geodesic is a straight line)*)
	\[Phi]y = PDTranslation[-x, y];
	
	(*Linear interpolation*)
	line = PDTranslation[x, # \[Phi]y] &/@ Subdivide[0,1, pts] // PDToXY;
	
	Return[ListLinePlot[line, PlotStyle->style,PlotRange->Full]];

);


PDInitCanvas[];

style = {Black}

polygon = {z0};
Iterate[]:= (
	polygon = polygon ~ Join ~ { Moebius[g["A"],Last[polygon]]};
);
Iterate[]&/@ Range[p-1];
polygonpts = PDHyperbolicPoint[#]&/@polygon;
polygongeos = Table[ PDGeodesic[polygon[[i]],polygon[[i+1]],50] , {i,1,7}]~ Join ~{PDGeodesic[polygon[[8]],polygon[[1]],20]};


fd = {0, z0, GeodesicMidpoint[z0, Moebius[g["A"],z0]] };
ptsfd = PDHyperbolicPoint[#]&/@fd;
geofd = { PDGeo[fd[[1]], fd[[2]],30,style],PDGeo[fd[[2]], fd[[3]],30,style],PDGeo[fd[[3]], fd[[1]],30,style] };

PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo;
PDUpdateCanvas[#]&/@geofd;
PDUpdateCanvas[#]&/@ptsfd;
gfx=PDShow[]


n = 8;
polygonT = Moebius[\[Gamma]B[#],z0] &/@ Range[1,8] //N ;
geosT = Table[ PDGeodesic[z0, polygonT[[i]],100] , {i,1,8}] ;
ptsT = PDHyperbolicPoint[#]&/@polygonT;
PDInitCanvas[];
(*PDUpdateCanvas[#]&/@geos;
PDUpdateCanvas[#]&/@pts;*)
PDUpdateCanvas[#]&/@ptsT
PDUpdateCanvas[#]&/@geosT;
gfx=PDShow[]


\[Gamma]1 = g["A"] . g["B"] . MatrixPower[g["A"],p/2 ] ;
\[Gamma][n_]:= MatrixPower[g["A"], n-1 ] . \[Gamma]1 . MatrixPower[g["A"], -n+1 ];


polygon = {z0};
Iterate[]:= (
	polygon = polygon ~ Join ~ { Moebius[g["A"],Last[polygon]]};
);
Iterate[]&/@ Range[p-1];
polygonpts = PDHyperbolicPoint[#]&/@polygon;


polygon


Moebius[\[Gamma][1],#]&/@ polygon


Translates[poly_]:=Table[Moebius[\[Gamma][i],#]&/@ poly,{i,1,8}] //N


polys ={polygon};



distance[z1_,z2_]:=ArcCosh[1 + (2 Abs[z1-z2]^2)/((1-Abs[z1]^2)(1-Abs[z2]^2))]//N


polys =DeleteDuplicates[polys~Join ~ArrayFlatten[Translates[#]&/@ polys,1]];
polyz = DeleteDuplicates[Flatten[polys]];
Length[polys]


polypts = PDHyperbolicPoint[#]&/@polyz;
Length[polypts]


DeleteDuplicates[polypts]//Length


geodata=SparseArray[DeleteCases[Flatten[Table[ If[ distance[polyz[[i]], polyz[[j]]]<3.1 && i!= j , {i,j}->1, None],{i,1, Length[polypts]},{j,1,i}]],None]]


polygeo = PDGeodesic[polyz[[#[[1]]]],polyz[[#[[2]]]],20] &/@  geodata["ExplicitPositions"];


PDInitCanvas[];
(*PDUpdateCanvas[#]&/@geos;
PDUpdateCanvas[#]&/@pts;*)

PDUpdateCanvas[#]&/@polypts;
PDUpdateCanvas[#]&/@polygeo;
gfx=PDShow[]


a[n_]:= Inverse[ \[Gamma][4n -1] . MatrixPower[g["A"],p/2 +2]]
b[n_]:= \[Gamma][4n ] . MatrixPower[g["A"],p/2 +2]


c[n_] := a[n] . b[n] . Inverse[a[n]] . Inverse[b[n]]


c[1] . c[2]//N // MatrixForm
