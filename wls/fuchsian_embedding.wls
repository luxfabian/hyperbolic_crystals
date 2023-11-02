#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Text:: *)
(*This script is only for the {8,3} triangle group*)


SetDirectory[NotebookDirectory[]];
 
Import["poincare.wls"];


p=8;
q=3;
n=-8;

If[n>0,
fprefix = "../data/"<>ToString[p]<>"_"<>ToString[q]<>"_modulo_"<>ToString[n],
fprefix = "../data/"<>ToString[p]<>"_"<>ToString[q]<>"_open_"<>ToString[-n]];

{dim}=Import[fprefix<>".info", "Data"]//ToExpression;

Print["Dimension of Hilbert space: ", dim]


(* ::Text:: *)
(*Read the basis file*)


lexikon = Import[fprefix<>".words", "Data"];
lexikon =#[[2]]&/@ SortBy[lexikon,First];


Adata = Import[fprefix<>"_A.reg","Table"];
Bdata = Import[fprefix<>"_B.reg","Table"];
iAdata = Import[fprefix<>"_iA.reg","Table"];
iBdata = Import[fprefix<>"_iB.reg","Table"];

\[CapitalPi]r["A"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Adata, Null], {dim,dim} ] ;
\[CapitalPi]r["B"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Bdata, Null], {dim,dim} ] ;
\[CapitalPi]r["iA"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iAdata, Null], {dim,dim} ] ;
\[CapitalPi]r["iB"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iBdata, Null], {dim,dim} ] ;


Adata = Import[fprefix<>"_A.lreg","Table"];
Bdata = Import[fprefix<>"_B.lreg","Table"];
iAdata = Import[fprefix<>"_iA.lreg","Table"];
iBdata = Import[fprefix<>"_iB.lreg","Table"];

\[CapitalPi]l["iA"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Adata, Null], {dim,dim} ] ;
\[CapitalPi]l["iB"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Bdata, Null], {dim,dim} ] ;
\[CapitalPi]l["A"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iAdata, Null], {dim,dim} ] ;
\[CapitalPi]l["B"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iBdata, Null], {dim,dim} ] ;


(* ::Text:: *)
(*Setup the generators*)


\[Alpha] =2 \[Pi]/p;
\[Beta]=2 \[Pi]/q;
r0  = Sqrt[Cos[(\[Alpha]+\[Beta])/2]/Cos[(\[Alpha]-\[Beta])/2]];

g["A"] =  { { Exp[I \[Alpha] /2],0},{0, Exp[-I \[Alpha]/2]}} ;
g["B"] = 1/(1-r0^2) { { Exp[I \[Beta] /2]-r0^2 Exp[-I \[Beta]/2],r0 (1-Exp[I \[Beta]]) Exp[I (\[Alpha]-\[Beta])/2]},{r0 (1-Exp[-I \[Beta]]) Exp[-I (\[Alpha]-\[Beta])/2],Exp[-I \[Beta]/2]-r0^2 Exp[I \[Beta]/2]}} ;


Moebius[ {{a_,b_},{c_,d_}}, z_]:= (a z + b )/(c z +d )


ParseWord[word_]:=(
z = z0;
( z = Moebius[g[#],z])&/@ Reverse[ Characters[ word ] ];
Return[z]
)


z0 =  Exp[I (\[Alpha]/2 + 0 \[Pi]/6)](r0 + 0.12 Exp[I \[Pi]]);

zs= N[ParseWord[#]]&/@lexikon;
pts = PDHyperbolicPoint[#]&/@zs;

z0 =  Exp[I (\[Alpha]/2 )] r0;
zsregular= N[ParseWord[#]]&/@lexikon;
ptsregular = PDHyperbolicPoint[#]&/@zsregular;


PDInitCanvas[];
PDUpdateCanvas[#]&/@pts;
PDUpdateCanvas[#]&/@{PDHyperbolicPoint[z0]};
gfx=PDShow[]


Ageo = PDGeodesic[zs[[#[[1]]]],zs[[#[[2]]]],20] &/@  \[CapitalPi]r["A"]["ExplicitPositions"];
Bgeo = PDGeodesic[zs[[#[[1]]]],zs[[#[[2]]]],20] &/@  \[CapitalPi]r["B"]["ExplicitPositions"];


AgeoRegular = PDGeodesic[zsregular[[#[[1]]]],zsregular[[#[[2]]]],20, {Thick,Orange}] &/@  \[CapitalPi]r["A"]["ExplicitPositions"];
BgeoRegular = PDGeodesic[zsregular[[#[[1]]]],zsregular[[#[[2]]]],20, {Thick,Orange}] &/@  \[CapitalPi]r["B"]["ExplicitPositions"];


PDInitCanvas[];
PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo;
PDUpdateCanvas[#]&/@AgeoRegular;
PDUpdateCanvas[#]&/@BgeoRegular; 
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


\[Alpha]f =2 \[Pi]/8;
\[Beta]f=2 \[Pi]/8;
r0f  = Sqrt[Cos[(\[Alpha]f+\[Beta]f)/2]/Cos[(\[Alpha]f-\[Beta]f)/2]];
fuchsian= Exp[I \[Alpha]f #] Exp[I \[Alpha]f/2]r0f &/@ Range[8]//N ;
fuchscolor=RGBColor["#9f0a0a"];
fuchsiangeo=PDGeodesic[fuchsian[[#]],fuchsian[[ Mod[#,8]+1]],50, {Thick,Black}]&/@ Range[8];


PDInitCanvas[];
PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo;
PDUpdateCanvas[#]&/@AgeoRegular;
PDUpdateCanvas[#]&/@BgeoRegular; 
PDUpdateCanvas[#]&/@fuchsiangeo; 
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


\[Gamma]1 = \[CapitalPi]l[#]&/@ { "A", "B","A","A","A","A","A", "B","A","A","A","A"};


\[Gamma]word = "ABAAAAABAAAA";


A[n_]:= StringJoin@@ ("A" &/@ Range[n])
iA[n_]:= StringJoin@@ ("D" &/@ Range[n])


\[CapitalPi]r["D"]=\[CapitalPi]r["iA"];


StringReplace[]


\[Gamma]w[n_]:=StringReplace[A[n-1]<>\[Gamma]word<>iA[n-1]]


\[Gamma]w[n_]:=StringReplace[A[n-1]<>\[Gamma]word<>iA[n-1],"AAAAAAAA"->"s"]


\[Gamma]w[n_]:=A[n-1]<>\[Gamma]word<>iA[n-1]


entireword = StringReplace[\[Gamma]w[1]<>\[Gamma]w[2+4]<>\[Gamma]w[3]<>\[Gamma]w[4+4]<>\[Gamma]w[1+4]<>\[Gamma]w[2]<>\[Gamma]w[3+4]<>\[Gamma]w[4],"AAAAAAAA"->""]


entireword =  StringReplace[StringReplace[StringReplace[StringReplace[\[Gamma]w[1]<>\[Gamma]w[2+4]<>\[Gamma]w[3]<>\[Gamma]w[4+4]<>\[Gamma]w[1+4]<>\[Gamma]w[2]<>\[Gamma]w[3+4]<>\[Gamma]w[4], {"DA"-> "","AD"-> ""}], {"DA"-> "","AD"-> ""}], {"DA"-> "","AD"-> ""}], {"DA"-> "","AD"-> ""}]


StringCount[ StringReplace[ entireword ,"AAAAAAAA"->"s"],"s"]


StringCount[entireword,"s"]


entireword/.{ "AAAA" -}


entireword = StringJoin@@ ("BAAAAABAA"&/@ Range[8])


loop = \[CapitalPi]r[#]&/@Characters[entireword] ;


seedword = 8;
Print["Seed word: ", lexikon[[ seedword ]] ];

vec = Table[ KroneckerDelta[i,seedword], {i,1,dim}];

Clear[\[Gamma]path];
\[Gamma]path[0]:= vec;
\[Gamma]path[i_]:= \[Gamma]path[i]=  Reverse[loop][[i]]  . \[Gamma]path[i-1];

orbit = \[Gamma]path[#] . zs &/@ Range[0,Length[loop]];
geoorbit =PDGeodesic[orbit[[#]],orbit[[ #+1]],300, {Thickness[0.008],fuchscolor}]&/@ Range[Length[orbit]-1];

\[Gamma]geo = PDGeodesicArrow[0 ,Moebius[ \[Gamma][#],0] ,300, {Thick,Gray}] &/@Range[8];
PDInitCanvas[];
PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo;
PDUpdateCanvas[#]&/@AgeoRegular;
PDUpdateCanvas[#]&/@BgeoRegular; 
PDUpdateCanvas[#]&/@fuchsiangeo; 
PDUpdateCanvas[#]&/@geoorbit; 
PDUpdateCanvas[#]&/@pts;
PDUpdateCanvas[#]&/@\[Gamma]geo; 
gfx=PDShow[]


Export["8_3_lattice.png",gfx, ImageResolution -> 300, Background -> None]


StringReplace[\[Gamma]w[2],"AAAAAAAA"->""]





\[Gamma]moeb=MatrixPower[ g["A"] . g["B"] . MatrixPower[g["A"],4], 2]//N;
\[Gamma][n_]:=MatrixPower[ g["A"], n-1] . \[Gamma]moeb . MatrixPower[ g["A"], -n+1]


\[Gamma]geo = PDGeodesicArrow[0 ,Moebius[ \[Gamma][#],0] ,300, {Thick,Gray}] &/@Range[8];
PDInitCanvas[];

PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo;
PDUpdateCanvas[#]&/@AgeoRegular;
PDUpdateCanvas[#]&/@BgeoRegular; 
PDUpdateCanvas[#]&/@fuchsiangeo; 
PDUpdateCanvas[#]&/@\[Gamma]geo; 
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


seedword = 5;
Print["Seed word: ", lexikon[[ seedword ]] ];

vec = Table[ KroneckerDelta[i,seedword], {i,1,dim}];

Clear[\[Gamma]path];
\[Gamma]path[0]:= vec;
\[Gamma]path[i_]:= \[Gamma]path[i]=  Reverse[\[Gamma]1][[i]]  . \[Gamma]path[i-1];

orbit = \[Gamma]path[#] . zs &/@ Range[0,Length[\[Gamma]1]];
geoorbit =PDGeodesic[orbit[[#]],orbit[[ #+1]],300, {Thick,fuchscolor}]&/@ Range[Length[orbit]-1];

\[Gamma]geo = PDGeodesicArrow[0 ,Moebius[ \[Gamma][#],0] ,300, {Thick,Gray}] &/@Range[8];
PDInitCanvas[];
PDUpdateCanvas[#]&/@Ageo;
PDUpdateCanvas[#]&/@Bgeo;
PDUpdateCanvas[#]&/@AgeoRegular;
PDUpdateCanvas[#]&/@BgeoRegular; 
PDUpdateCanvas[#]&/@fuchsiangeo; 
PDUpdateCanvas[#]&/@\[Gamma]geo; 
PDUpdateCanvas[#]&/@geoorbit; 
PDUpdateCanvas[#]&/@pts;
gfx=PDShow[]


Clear[\[Alpha],\[Beta]]


Aa = 8 \[Pi] - 8 \[Alpha] - 2 \[Pi] ;
Ab = 3 \[Pi] - 3 \[Beta] - 2 \[Pi] ;
Ac = - 2 (\[Pi]-\[Alpha]-\[Beta]);


A\[Gamma] = 8 ( Aa/8 + 2 Ab + 4 Ac)


Agauge = - (\[Pi] - \[Alpha]);


Bgauge = - (\[Pi] - \[Beta]);


3 Aa -A\[Gamma] + 2 8 Bgauge + 7 8 Agauge  //FullSimplify // Expand // TeXForm


3 Aa - A\[Gamma] //FullSimplify
