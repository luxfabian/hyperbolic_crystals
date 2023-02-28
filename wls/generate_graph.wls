#!/usr/bin/env wolframscript
(* ::Package:: *)

dir = If[$Notebooks, NotebookDirectory[], Directory[] ];
SetDirectory[dir];


p=5;
q=4;


fname = ToString[p]<>"_"<>ToString[q]<>"_A.reg"


Aname = "../data/5_4_modulo_2_A.reg";
Bname = "../data/5_4_modulo_2_B.reg";
ABname = "../data/5_4_modulo_2_AB.reg";


Aname = "../data/5_4_open_5_A.reg";
Bname = "../data/5_4_open_5_B.reg";
ABname = "../data/5_4_open_5_AB.reg";


Adata = Import[Aname,"Table"];
Bdata = Import[Bname,"Table"];
ABdata = Import[ABname,"Table"];


Style[ DirectedEdge[i , ctab[[i]]],colors[op]]


Agraph = Flatten[ Style[DirectedEdge[ #[[1]] ,#[[2]]],  Blue] &/@ Adata ];
Bgraph = Flatten[Style[DirectedEdge[ #[[1]] ,#[[2]]],  Red] &/@ Bdata ];
Cgraph = Flatten[Style[DirectedEdge[ #[[1]] ,#[[2]]],  Green] &/@ ABdata ];


graph = Agraph ~ Join ~ Bgraph ~ Join ~ Cgraph;


GraphPlot3D[graph]


Jopi


graph =Graph[Flatten[ Agraph]];


Agraph[[4]]


GraphPlot[Flatten[Bgraph]]


GraphPlot3D[]
