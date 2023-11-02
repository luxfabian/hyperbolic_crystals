#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[NotebookDirectory[]];

Import["poincare.wls"];


PDInitCanvas[];

id = 0;
g1 = 0.75I;
g2 = Exp[2 \[Pi] I /3] g1;
g3 = Exp[2 \[Pi] I /3] g2;

geostyle = {Black};
geopts = 50;
geos = { PDGeodesic[id,g1,geopts,geostyle],PDGeodesic[id,g2,geopts,geostyle],PDGeodesic[id,g3,geopts,geostyle],PDGeodesic[g1,g2,geopts,geostyle]
,PDGeodesic[g2,g3,geopts,geostyle],PDGeodesic[g3,g1,geopts,geostyle]};
pts = PDHyperbolicPoint[#]&/@{id,g1,g2,g3};


PDUpdateCanvas[#]&/@ geos;
PDUpdateCanvas[#]&/@ pts;

gfx=PDShow[]
