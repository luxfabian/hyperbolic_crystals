#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[If[$Notebooks, NotebookDirectory[], Directory[] ]];


data=Import["./gapout_fuchsian/coset_tables.csv", "Text"];
CosetTable =ToExpression[ StringReplace[data,{ "[" -> "{", "]" -> "}"} ] ] ;


data = Import["./gapout_14/quotients.csv","Text"];
Classification =ToExpression[ StringReplace[data,{ "[" -> "{", "]" -> "}"} ] ] ;


Length[Classification]


Length[#[[1]]] &/@CosetTable


coherent = {{5,160},{7,360}, {15,2560}, {42, 15000 }};


lp1 = ListPlot[ Length[#[[1]]] &/@CosetTable, GridLines->Automatic, Frame-> True, FrameLabel->{"normal subgroup","index"} ];
lp2 = ListPlot[ coherent, PlotStyle-> Red,PlotMarkers->{Automatic, 6 }];
Show[lp1, lp2,ImageSize-> 400]


bc=BarChart[
	Sort[Counts[Classification]],  
	ChartLabels -> Automatic, 
	BarOrigin->Left,
	LabelStyle->12,
	ImageSize->500, 
	AspectRatio->2, 
	ColorFunction->"DarkBands", 
	GridLines->Automatic,
	AxesLabel->"Abs. frequency"
]
Export["./gapout/taxonomy.png", bc, ImageResolution -> 300]


(* ::Text:: *)
(*Index of group generators (position of generators in basis set)*)


a = 1; 
b = 3;


colors[a]:= Red;
colors[b]:= Blue;
Represent[op_, sg_]:=(
	ctab = CosetTable[[sg]][[op]];
	graph = Table[Style[ DirectedEdge[i , ctab[[i]]],colors[op]], {i,1, Length[ctab]} ];
	Return[graph];
);


PureGraph[sg_]:=Represent[a,sg]~ Join ~ Represent[b,sg]
Cayley[sg_]:= GraphPlot[PureGraph[sg] ,GraphLayout->"SpringElectricalEmbedding"]
Cayley3D[sg_]:= GraphPlot3D[ PureGraph[sg],GraphLayout->"SpringElectricalEmbedding"]


Represent[b,4]//GraphPlot


(* ::Text:: *)
(*Matrix representation of the generators*)


MatRepresent[op_, sg_]:=(
	ctab = CosetTable[[sg]][[op]];
	
	dim = Length[ctab]; 
	
	mat = Table[ KroneckerDelta[i, ctab[[j]]],  {i,1,dim},{j,1,dim}];
	
	Return[mat];
);


(* ::Text:: *)
(*Definition of the Hamiltonian:*)


Hamiltonian[sg_]:= Sum[ MatRepresent[i,sg],{i,1,4}]/4


(* ::Text:: *)
(*and its spectrum:*)


Spectrum[sg_] :=Sort[Eigenvalues[N[Hamiltonian[sg]]]]


(* ::Text:: *)
(*The integrated density of states can be computed from the spectrum as*)


IDS[\[Mu]_, spec_]:= N[ Count[  # < \[Mu] &/@ spec, True]/ Length[spec] ] 


(* ::Text:: *)
(*For the purposes of plotting the IDS, we define*)


pts = 1000;
\[Mu]s = Subdivide[-1,1,pts]//N;

ids[spec_]:= Transpose[{\[Mu]s, IDS[#,spec]&/@\[Mu]s}];


PlotSpectrum[sg_]:=ListLinePlot[
	ids[Spectrum[sg]] ,
	Frame-> True, 
	FrameLabel->{"x","F(x)"},
	(*PlotLabel-> "\!\(\*SubscriptBox[\(F\), \(2\)]\)/H = "<>Classification[[sg]],*)
	FrameStyle->20,
	LabelStyle->20
] 


PlotSpectrum[5]


Analyze[sg_]:=(
	gfx = GraphicsGrid[{{ Cayley[sg]}, {PlotSpectrum[sg] }}];
	Export["./gapout_14/"<>StringPadLeft[ToString[sg],4,"0"]<>".png",gfx,ImageResolution-> 300];
	Return[gfx];
)


Cayley[6]


gfx=Analyze[#]&/@ Range[Length[CosetTable]];


(* ::Text:: *)
(*Indices for which the factor group is isomorphic to SL(2,3)*)


SL23indices= DeleteCases[If[  Classification[[#]]== "SL(2,3)", #,None]&/@ Range[Length[CosetTable]] ,None]


(* ::Text:: *)
(*Adapt plotting for this case:*)


PlotSpectrum[sg_]:=ListLinePlot[
	ids[Spectrum[sg]] ,
	Frame-> True, 
	FrameLabel->{"x","F(x)"},
	FrameStyle->30,
	LabelStyle->30,
	PlotLabel-> Style[ Classification[[sg]], 22]
] 

Analyze[sg_]:=(
	gfx = GraphicsGrid[{{ CayleyGraphF2[sg]}, {PlotSpectrum[sg] }}];
	Return[gfx];
)


(* ::Text:: *)
(*Choose random subset*)


indices = RandomInteger[{1,Length[CosetTable]}] &/@ Range[16] 


(* ::Text:: *)
(*Except for the following case (gives reproducible figure for index 32)*)


If[ Length[CosetTable] == 151,
indices = {58,11,45,14,24,84,110,47,59,23,144,28,142,86,102,15}];


gfx=Analyze[#]&/@ indices;


spec=Grid[ ArrayReshape[gfx,{4,4}] ];


Export["./gapout/kaleidoscope.png",spec,ImageResolution-> 300];
