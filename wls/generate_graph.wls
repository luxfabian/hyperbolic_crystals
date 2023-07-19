#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*generate_graph.wls*)


(* ::Text:: *)
(*Author: Fabian R. Lux*)
(*Date: 7/19/23*)
(**)
(*Reads the regular representation files and constructs the Cayley graph. Can be used to verify the validity of the generated datasets*)


dir = If[$Notebooks, NotebookDirectory[], Directory[] ];
SetDirectory[dir];


Aname = "../data/5_4_modulo_2_A.reg";
Bname = "../data/5_4_modulo_2_B.reg";
ABname = "../data/5_4_modulo_2_AB.reg";

Adata = Import[Aname,"Table"];
Bdata = Import[Bname,"Table"];
ABdata = Import[ABname,"Table"];

Acolor = RGBColor["#4472c4"]
Bcolor = RGBColor["#ffc000"]
Ccolor = RGBColor["#d42525"]


Agraph = Flatten[ Style[DirectedEdge[ #[[1]] ,#[[2]]],  Acolor] &/@ Adata ];
Bgraph = Flatten[Style[DirectedEdge[ #[[1]] ,#[[2]]],  Bcolor] &/@ Bdata ];
Cgraph = DeleteCases[Flatten[ If[ #[[1]] <#[[2]],Style[UndirectedEdge[ #[[1]] ,#[[2]]],  Ccolor]] &/@ ABdata ],Null];
graph = Graph[Agraph ~ Join ~ Bgraph ];


GraphPlot[Agraph,VertexStyle->Black]


GraphPlot[Bgraph,VertexStyle->Black]


GraphPlot[Cgraph,VertexStyle->Black]


graph =  Graph[ Agraph ~ Join ~ Bgraph ~ Join ~ Cgraph];


opts = GraphLayout->{"SpringElectricalEmbedding", "SpringConstant"-> 1,"RepulsiveForcePower"->-1, "Tolerance"->10^-18,"MaxIteration"->10^6};


GraphPlot[graph, VertexStyle->Black, opts]
