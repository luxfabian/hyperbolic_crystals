#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[NotebookDirectory[]];

Import["poincare.wls"];


(* ::Text:: *)
(*Get info about group*)


{g,n} = #[[2]]&/@Import["../data/tffg_group_specs.inp", "Table"];

If[n>0,
fprefix = "../data/tffg_"<>ToString[g]<>"_modulo_"<>ToString[n],
fprefix = "../data/tffg_"<>ToString[g]<>"_open_"<>ToString[-n]];

{dim}=Import[fprefix<>".info", "Data"]//ToExpression;

Print["Dimension of Hilbert space: ", dim]


(* ::Text:: *)
(*Load right-regular representation of generators*)


Clear[reg,arr,graph];
reg[i_]:=reg[i]=Import[fprefix<>"_"<>ToString[i-1]<>".reg","Table"];
arr[i_]:= arr[i]=SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ reg[i], Null], {dim,dim} ] ;
graph[i_]:=  graph[i] = DirectedEdge[#[[1]],#[[2]]] &/@  arr[i]["ExplicitPositions"];


graphs = graph[#]&/@ Range[2g];


MatrixPower[arr[4],12] == IdentityMatrix[dim]


GraphPlot[graphs[[2]]]


GraphPlot[ Join@@ graphs[[1;;2g ]] ]


H = Sum[arr[i],{i,1, 4 g} ];
spec0 = Eigenvalues[N[Normal[H]]];


ListPlot[Sort[spec0]]


\[Mu]s = Subdivide[-6,6,500];

ids[\[Mu]_]:=1/Length[spec0] Plus@@ ( If[# < \[Mu] ,1,0] &/@ spec0  ) //N

ids0 = ids[#] &/@ \[Mu]s ;

ListLinePlot[Transpose[{\[Mu]s, ids0}]]


spec1 = #[[2]] &/@ magspec[1];
ListPlot[Sort[spec1]]


spec2 = #[[2]] &/@ magspec[2];
ListPlot[Sort[spec2]]


spec3 = #[[2]] &/@ magspec[3];
ListPlot[Sort[spec3]]





4 dim


ListPlot[{Sort[spec1],Sort[spec2]}]


Mod[i+1,2]+1


magspec[k_]:=(
	s = Exp[(2 \[Pi] I)/k];
	
	Clear[\[CapitalSigma]];
	\[CapitalSigma][i_/;i>2]:= SparseArray[ IdentityMatrix[k] ];
	\[CapitalSigma][1] := If[k>1,SparseArray[{Band[{2,1}]->1},{k,k}], SparseArray[{},{k,k}]]  + SparseArray[{ {1,k}-> 1}, {k,k}];
	\[CapitalSigma][2] := SparseArray[ {#,#}->s^(#-1) &/@ Range[k], {k,k}];
	\[CapitalSigma][1 + 2 g] := Inverse[\[CapitalSigma][1]];
	\[CapitalSigma][2 + 2 g] := Inverse[\[CapitalSigma][2]];
	
	Hmag = Sum[ KroneckerProduct[ \[CapitalSigma][i] , arr[i]],{i,1,4g } ] // N ;
	Hmagdense = Normal[Hmag];
	
	spec = {k/dim,#}&/@ Eigenvalues[N[Hmagdense]];
	
	Return[spec];
);


KroneckerProduct[ \[CapitalSigma][i] , arr[i]]


\[Gamma][i_]:= arr[i]

k=4;
s = Exp[(2 \[Pi] I)/k];

Clear[\[CapitalSigma]];
\[CapitalSigma][i_/;i>2]:= SparseArray[ IdentityMatrix[k] ];
\[CapitalSigma][1] := If[k>1,SparseArray[{Band[{2,1}]->1},{k,k}], SparseArray[{},{k,k}]]  + SparseArray[{ {1,k}-> 1}, {k,k}];
\[CapitalSigma][2] := SparseArray[ {#,#}->s^(#-1) &/@ Range[k], {k,k}];
\[CapitalSigma][1 + 2 g] := Inverse[\[CapitalSigma][1]];
\[CapitalSigma][2 + 2 g] := Inverse[\[CapitalSigma][2]];

\[Gamma]mag[i_]:= KroneckerProduct[ \[CapitalSigma][i],arr[i]];


MatrixForm[\[CapitalSigma][1]]


MatrixForm[\[CapitalSigma][2]]


MatrixForm[\[CapitalSigma][3]]


\[Gamma][1] . \[Gamma][2] . \[Gamma][3] . \[Gamma][4] . \[Gamma][5] . \[Gamma][6] . \[Gamma][7] . \[Gamma][8] == IdentityMatrix[dim]


(s \[Gamma]mag[1] . \[Gamma]mag[2] . \[Gamma]mag[3] . \[Gamma]mag[4] . \[Gamma]mag[5] . \[Gamma]mag[6] . \[Gamma]mag[7] . \[Gamma]mag[8] )== IdentityMatrix[dim k]


op =  \[Gamma][1] . \[Gamma][2] . \[Gamma][3] . \[Gamma][4] . \[Gamma][5] . \[Gamma][6] . \[Gamma][7] . \[Gamma][8] ;


MatrixPower[ op,k]== -IdentityMatrix[dim k]


\[Gamma][1] . \[Gamma][2] . \[Gamma][3] . \[Gamma][4] . \[Gamma][5] . \[Gamma][6] . \[Gamma][7] . \[Gamma][8] == IdentityMatrix[dim]


spec4 = #[[2]] &/@ magspec[3];
ListPlot[Sort[spec4]]


Mod[1+1,2]+1


k = 3;

s = Exp[(2 \[Pi] I)/k];
	
Clear[\[CapitalSigma]];
\[CapitalSigma][i_]:= SparseArray[ IdentityMatrix[k] ];
\[CapitalSigma][1] := If[k>1,SparseArray[{Band[{2,1}]->1},{k,k}], SparseArray[{},{k,k}]]  + SparseArray[{ {1,k}-> 1}, {k,k}];
\[CapitalSigma][2] := SparseArray[ {#,#}->s^(k+#-1) &/@ Range[k], {k,k}]


s (\[CapitalSigma][1] . \[CapitalSigma][2] . Inverse[\[CapitalSigma][1]] . Inverse[\[CapitalSigma][2]]) //MatrixForm


k 


s


\[CapitalSigma][1]//MatrixForm


\[CapitalSigma][2]//MatrixForm


s^(k-1)


s (\[CapitalSigma][1] . \[CapitalSigma][2] . Inverse[\[CapitalSigma][1]] . Inverse[\[CapitalSigma][2]]) //MatrixForm





spec2= magspec[3];


spec1s = #[[2]] &/@ spec1;


spec2s = #[[2]] &/@ spec2;


spec1= magspec[2];
spec1s = #[[2]] &/@ spec1;

ids[\[Mu]_]:=1/Length[spec1s] Plus@@ ( If[# < \[Mu] ,1,0] &/@ spec1s  ) //N
 
ids1 = ids[#] &/@ \[Mu]s ;


spec2= magspec[2];
spec2s = #[[2]] &/@ spec2;

ids[\[Mu]_]:=1/Length[spec2s] Plus@@ ( If[# < \[Mu] ,1,0] &/@ spec2s  ) //N
 
ids2 = ids[#] &/@ \[Mu]s ;

ListLinePlot[Transpose[{\[Mu]s, ids2}]]


ids[\[Mu]_]:=1/Length[spec2s] Plus@@ ( If[# < \[Mu] ,1,0] &/@ spec2s  ) //N
 
ids2 = ids[#] &/@ \[Mu]s ;

ListLinePlot[Transpose[{\[Mu]s, ids2}]]


dim 3


ListLinePlot[{Transpose[{\[Mu]s, ids0}],Transpose[{\[Mu]s, ids1}],Transpose[{\[Mu]s, ids2}]}, PlotLegends->{"No magnetic field","k=2","k=3"},Frame-> True, FrameLabel->{"\[Mu]","IDS"}]


ids1


ks = Range[2,30];


dim dim


specs  = magspec[#]&/@ks;


ListPlot[ ArrayFlatten[specs,1] ]


DensityHistogram[ArrayFlatten[specs,1],{Round[Length[ks]/2],50}]


ArrayFlatten[specs,1]  //ListPlot


k = 2;
s = Exp[(2 \[Pi] I)/k];

Clear[\[CapitalSigma]]
\[CapitalSigma][i_]:= SparseArray[ IdentityMatrix[k] ];
\[CapitalSigma][1] := SparseArray[{Band[{2,1}]->1},{k,k}];
\[CapitalSigma][2] := SparseArray[ {#,#}->s^(k+#-1) &/@ Range[k], {k,k}];

Hmag = Sum[ KroneckerProduct[ \[CapitalSigma][i] , arr[i]],{i,1,2g } ] // N ;
Hmag = Hmag + ConjugateTranspose[Hmag];
Hmagdense = Normal[Hmag];

spec = Eigenvalues[N[Hmagdense]];

\[Mu]s = Subdivide[-6,6,500];
idss = ids[#] &/@ \[Mu]s ;

ListLinePlot[Transpose[{\[Mu]s, idss}]]





Hdense = Normal[H] ;


spec = Eigenvalues[N[Hdense]];


spec = Eigenvalues[N[Hmagdense]];


ids[\[Mu]_]:=1/Length[spec] Plus@@ ( If[# < \[Mu] ,1,0] &/@ spec  ) //N


\[Mu]s = Subdivide[-6,6,500];
idss = ids[#] &/@ \[Mu]s ;


ListLinePlot[Transpose[{\[Mu]s, idss}]]


ListPlot[Sort[spec]]


Eigenvalues[H]


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


Ahlfors[z1_,z2_] := Sqrt[(1 - Abs[z1]^2)(1 - Abs[z2]^2) + Abs[z1-z2]^2]


GeodesicMidpoint[z1_, z2_]:=(z2 (1- Abs[z1]^2)+z1 (1- Abs[z2]^2))/(1- Abs[z1]^2  Abs[z2]^2+Ahlfors[z1,z2] Sqrt[(1-Abs[z1]^2)(1-Abs[z2]^2)])


(* ::Text:: *)
(*Seed coordinate*)


z0 =  Exp[I \[Alpha]/2]r0;


(* ::Text:: *)
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
geofd = { PDGeodesic[fd[[1]], fd[[2]],30,style],PDGeodesic[fd[[2]], fd[[3]],30,style],PDGeodesic[fd[[3]], fd[[1]],30,style] };

PDUpdateCanvas[#]&/@geofd;

PDUpdateCanvas[#]&/@ptsdual;

gfx=PDShow[]





genus = p/4 


(* ::Text:: *)
(*Old definition:*)


\[Gamma][1] := g["A"] . g["B"] . MatrixPower[g["A"],2 genus]
\[Gamma][n_]:= MatrixPower[g["A"], (n-1)] . \[Gamma][1] . MatrixPower[g["A"], -(n-1)]
\[Gamma]s = Table[\[Gamma][i],{i,1, }] ;
Dot@@ \[Gamma]s  //N //Chop//MatrixForm


(* ::Text:: *)
(*New definition :*)


\[Gamma][1] := g["A"] . g["B"] . MatrixPower[g["A"],2 genus]
\[Gamma][n_]:= MatrixPower[g["A"],-(2 genus-1 ) (n-1)] . \[Gamma][1] . MatrixPower[g["A"],(2 genus-1 ) (n-1)]
\[Gamma]s = Table[\[Gamma][i],{i,1,p/2}] ;
Dot@@ \[Gamma]s  //N //Chop//MatrixForm








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


magspec[k_, n_]:=(
	s = Exp[n (2 \[Pi] I)/k];
	
	Clear[\[CapitalSigma]];
	\[CapitalSigma][i_/;i>2]:= SparseArray[ IdentityMatrix[k] ];
	\[CapitalSigma][1] := If[k>1,SparseArray[{Band[{2,1}]->1},{k,k}], SparseArray[{},{k,k}]]  + SparseArray[{ {1,k}-> 1}, {k,k}];
	\[CapitalSigma][2] := SparseArray[ {#,#}->s^(#-1) &/@ Range[k], {k,k}];
	\[CapitalSigma][1 + 2 g] := Inverse[\[CapitalSigma][1]];
	\[CapitalSigma][2 + 2 g] := Inverse[\[CapitalSigma][2]];
	
	Hmag = Sum[ KroneckerProduct[ \[CapitalSigma][i] , arr[i]],{i,1,4g } ] // N ;
	Hmagdense = Normal[Hmag];
	
	spec = {n/k,#}&/@ Eigenvalues[N[Hmagdense]];
	
	Return[spec];
);


kv = 40;
ns = Range[0,kv];
specs = magspec[kv,#] &/@ ns;


ListPlot[ ArrayFlatten[specs,1], Frame-> True, GridLines->Automatic, FrameLabel->{" Flux \[Phi] ( s =  \!\(\*SuperscriptBox[\(e\), \(\(2\)  \(\[Pi]\)\(\\\ \)\(\[ImaginaryI]\)\(\\\ \)\(\[Phi]\)\(\\\ \)\)]\))", "Energy"},AspectRatio->1]
