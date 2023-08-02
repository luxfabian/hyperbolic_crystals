#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Poincar\[EAcute] disk model of the hyperbolic plane*)


(* ::Subtitle:: *)
(*Collection of common operations*)


(* ::Text:: *)
(*Author: 	Fabian R. Lux*)
(*Date:	1/30/2023*)
(**)
(*This Wolfram Language Script implements geometric manipulations of points on the Poincar\[EAcute]	disk model of the hyperbolic plane*)


(* ::Chapter:: *)
(*Canvas and color schemes*)


(* ::Text:: *)
(*The default color scheme is defined by*)


PDVertexColor =RGBColor[  "#1e81b0" ];
PDOuterColorOfPoincareDisk =  LightBlue;
PDInnerColorOfPoincareDisk =White;
PDGeodesicColor = Black;


PDSetGeodesicColor[col_]:=(
	PDGeodesicColor = col;
);

PDSetGeodesicColor::usage = "Set the color of a geodesic| on the Poincar\[EAcute] disk." ;


PDSetVertexColor[col_]:=(
	PDVertexColor = col;
);

PDSetVertexColor::usage = "Set the color of a point (a vertex) on the Poincar\[EAcute] disk." ;


PDSetOuterColorOfPoincareDisk[col_]:=(
	PDOuterColorOfPoincareDisk = col;
);

PDSetVertexColor::usage = "The Poincar\[EAcute] disk is colored with a radial gradient filling. This routine sets the outer color." ;


PDSetInnerColorOfPoincareDisk[col_]:=(
	PDInnerColorOfPoincareDisk = col;
);

PDSetVertexColor::usage = "The Poincar\[EAcute] disk is colored with a radial gradient filling. This routine sets the inner color." ;


PDGetCanvas[]:=(
	PoincareDisk:=Graphics[{RadialGradientFilling[{PDInnerColorOfPoincareDisk,PDOuterColorOfPoincareDisk }],EdgeForm[Directive[Black]], Disk[{0,0},1],EdgeForm[Directive[Black]], Circle[{0,0},1]}];
	Return[{PoincareDisk}];
);

PDGetCanvas::usage = "Return a graphics object which represents the Poincar\[EAcute] disk.";


PDInitCanvas[]:=(
	PDcanvas = PDGetCanvas[];
)

PDInitCanvas::usage = "Defines the internal variable PDcanvas, which represents the graphics pipeline. It contains a representation of the Poincar\[EAcute] disk and all objects which should eventually be drawn onto it.";


PDUpdateCanvas[obj_]:=(
	PDcanvas = PDcanvas ~ Join ~ { obj };
);

PDUpdateCanvas::usage = "PDUpdateCanvas[obj] adds an object to the canvas.";


PDShow[]:=(
	Show[PDcanvas]
);

PDShow::usage = " PDShow[] will plot all the objects which have been added to the canvas.";


(* ::Chapter:: *)
(*Hyperbolic point*)


PDHyperbolicPoint[z_]:=(
	Clear[d];
	d = ArcTanh[Abs[z]];
	
	ptsize=0.03;
	
	radius = ptsize Exp[-d];
	point = Graphics[{EdgeForm[Thin],PDVertexColor, Disk[{Re[z],Im[z]},radius]}];
	
	Return[point];
);

PDHyperbolicPoint::usage = "Given a complex number z, PDHyperbolicPoint[z] will return a graphics object which injects the point into the hyperbolic plane.";


(* ::Chapter:: *)
(*M\[ODoubleDot]bius transformations*)


PDAutomorphism[\[Alpha]_, a_, z_]:= Exp[I \[Alpha]] (z-a)/(1- Conjugate[a] z)

PDAutomorphism::usage = "PDAutomorphism[\[Alpha], a, z] maps the complex number z to 0 for z=a and rotates the result by angle \[Alpha]. The point z=0 is mapped to -a.";


PDTranslation[a_, z_]:= PDAutomorphism[0, -a, z];

PDTranslation::usage =" PDTranslation[a,z] is a rigid translation in the hyperbolic plane which maps z=0 to a.";


PDRotation[\[Alpha]_, z_]:= PDAutomorphism[\[Alpha], 0, z]

PDRotation::usage="Rotate z by \[Alpha] around the point 0";


(* ::Chapter:: *)
(*Hyperbolic circle*)


PDCircle[pts_, a_, r_]:=(
	(*Create a circle around the origin*)
	
	rad = Sqrt[ ((Cosh[r]-1)/(Cosh[r]+1))] ; (*euclidean radius*)
	
	zs = PDTranslation[a, Exp[I 2 \[Pi] # / pts ] rad ] &/@ Range[pts+1];
	
	xydat = { Re[#] , Im[#] } &/@ zs;
	
	gfx = ListLinePlot[xydat, PlotStyle->{Thin,Black}];
	Return[gfx];
);

PDCircle::usage = "PDCircle[pts, a, r] creates a circle of hyperbolic radius r around the point a. pts controls the number of points used to interpolate the circle.";


(* ::Chapter:: *)
(*Geodesic*)


PDToXY[list_]:= {Re[#],Im[#]}&/@ list ;

PDToXY::usage = "Turns a list of complex values to a list of x and y Euclidean coordinates.";


PDGeodesic[x_, y_, pts_,style_:{Thin,PDGeodesicColor}]:=(

	(*Translate point x to 0 (then the geodesic is a straight line)*)
	\[Phi]y = PDTranslation[-x, y];
	
	(*Linear interpolation*)
	line = PDTranslation[x, # \[Phi]y] &/@ Subdivide[0,1, pts] // PDToXY;
	
	Return[ListLinePlot[line, PlotStyle->style,PlotRange->Full]];

);

PDGeodesic::usage = "Draw a geodesic between the complex points x and y. pts controls the linear interpolation.";


(* ::Chapter:: *)
(*Triangles*)


PDRightTrianglePoints[p_, q_, r_]:=(

	Q = 1/p+1/q+1/r<1;
	
	If[Not[Q], 
		{ 
			Print["\!\(\*
StyleBox[\"Impossible\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"Triangle\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\":\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"1\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"/\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"p\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"+\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"1\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"/\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"q\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"+\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"1\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"/\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"r\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\">\",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)\!\(\*
StyleBox[\"1\",\nFontWeight->\"Bold\"]\)"];
			Abort[];
		}
	];
	
	\[Alpha] := \[Pi]/p;
	\[Beta] := \[Pi]/q;
	\[Gamma] := \[Pi]/r;
	
	(* Law of cosines II *)
	coshc = ( Cos[\[Alpha]] Cos[\[Beta]] + Cos[\[Gamma]]) / ( Sin[\[Alpha]] Sin[\[Beta]] );
	c = ArcCosh[coshc];
	
	(*law of sines *)
	sinha = Sin[\[Alpha]] Sinh[c]/ Sin[\[Gamma]];
	a = ArcSinh[sinha];
	
	sinhb = Sin[\[Beta]] Sinh[c]/ Sin[\[Gamma]];
	b = ArcSinh[sinhb];
	
	
	rb =  Sqrt[ ((Cosh[b]-1)/(Cosh[b]+1))];
	
	r0 = Sqrt[(coshc-1)/(coshc+1)];
	ptA = 0;
	ptB = Exp[I \[Alpha]] r0;
	ptC = rb;
	
	Return[{ptA, ptB, ptC,a,b,c,\[Alpha],\[Beta],\[Gamma]}];
);


PDRightTrianglePoints::usage="PDRightTrianglePoints[p, q, r] return the complex coordinates of a triangle ABC in the hyperbolic plane with vertex A at 0 and vertex C on the positiv real axis. The angles at A,B,C are \[Pi]/p, \[Pi]/q and \[Pi]/r respevtively. One must have that 1/p+1/q+1/r < 1. pts controls the number of points used for the interpolation of geodesics.";


PDRightTriangle[p_, q_, r_, pts_]:=(

	{ptA, ptB, ptC,a,b,c,\[Alpha],\[Beta],\[Gamma]} = PDRightTrianglePoints[p,q,r];
	
	pipeline = {PDGeodesic[ptA,ptB, pts],PDGeodesic[ptB,ptC, pts],PDGeodesic[ptC,ptA, pts],PDHyperbolicPoint[ptA], PDHyperbolicPoint[ptB], PDHyperbolicPoint[ptC]} ;
	Return[pipeline];
);

PDRightTriangle::usage="PDRightTriangle[p, q, r, pts] generates a triangle ABC in the hyperbolic plane with vertex A at 0 and vertex C on the positiv real axis. The angles at A,B,C are \[Pi]/p, \[Pi]/q and \[Pi]/r respevtively. One must have that 1/p+1/q+1/r < 1. pts controls the number of points used for the interpolation of geodesics.";
