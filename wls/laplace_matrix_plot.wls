#!/usr/bin/env wolframscript
(* ::Package:: *)

dir = If[$Notebooks, NotebookDirectory[], Directory[] ];
SetDirectory[dir];


fix = "open_8";


fix = "modulo_4";


Aname = "../data/5_4_"<>fix<>"_A.reg";
Bname = "../data/5_4_"<>fix<>"_B.reg";
iAname = "../data/5_4_"<>fix<>"_iA.reg";
iBname = "../data/5_4_"<>fix<>"_iB.reg";
ABname = "../data/5_4_"<>fix<>"_AB.reg";
Specinfoname= "../data/5_4_"<>fix<>".info";


Adata = Import[Aname,"Table"];
Bdata = Import[Bname,"Table"];
iAdata = Import[iAname,"Table"];
iBdata = Import[iBname,"Table"];
ABdata = Import[ABname,"Table"];

d = ToExpression[Import[Specinfoname, "Data"][[1]]];

\[CapitalPi]["A"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Adata, Null], {d,d} ] ;
\[CapitalPi]["iA"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iAdata, Null], {d,d} ] ;
\[CapitalPi]["B"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ Bdata, Null], {d,d} ] ;
\[CapitalPi]["iB"]= SparseArray[ DeleteCases[ If[ #[[1]]>= 0 && #[[2]] >= 0 , ( {#[[1]]+1,#[[2]]+1} -> 1 ) ]  &/@ iBdata, Null], {d,d} ] ;


MatrixPlot[\[CapitalPi]["A"]+ \[CapitalPi]["iA"]+\[CapitalPi]["B"]+ \[CapitalPi]["iB"]]


MatrixPlot[\[CapitalPi]["A"]+ \[CapitalPi]["iA"]+\[CapitalPi]["B"]+ \[CapitalPi]["iB"]]


Graphics[Raster[CompressedData["
1:eJzt3X1wl9WVwHGyAlUglpeAinZxq2NhGCqYtVChudSxi9bMLKKDiG+8WBG7
spUXNXRqtCw1MYDIiLi04pIKDCUgbINFXQqMRWu7IgUqSujmDfL2+/2ChNIp
CHRnh870nz0Zjp7Lvc8v3850Pm3mmef3PPf5kefk3HPP/Ycp/zruOxd06tRp
6l//+3//++x/PnGIiIiISfXEvunjhozNuEE3rPzF0J82u9x5D/3m1furXe8l
B6/vty7tbvunHmt/8z8Nrnd13gc5Pfe5Q89sPDFxe8Z9e2fznX1vTbkXWlY+
Oe1MvXj8lKq1X/50+xH3p75LX/n+fVWuZnLD4xsu+3rw+0ZERMSO6+n8gV/8
aG/GvbW+ZtcbL6Tcl779L5P/uPWQO/Zaw6NHf5JxPcqvuOhbO1tc+fLSujvL
69y8zTl3z2/KuIt3D99dVZdylbu6jbu29LDrfcOigoXXtLquw7+x75vvp93s
exquqmlpdANH//GRG4ta3djt//7xe11a3a9nTGr47SUpVzRs+9YVmzLu/l/N
6jatT8p9fcDAvUseqXe521ZMGDs47TZXPJW7r6DWdZ6x5PD66U+I8VXo8UNE
RMTk+aVvfPPDxn9udeu+n/pw/vSMe+U7V3UZ/3SzK/9l/y8vfaDVnXn5xo/W
/CLj/rFlx/HFhS3ugcETFzXfl3H9Rl7SpfpAszvUY0PPNaNrxDhKOv9Lb017
9oaRrW7CXS/8furlGXdh7b/lXf56k1u/4/03jixsdT+t++RA6vZWd/+xMx92
rU+5C4dXrb7udMZNblv0u9ND0+7m29fe+PAdDa5tavn4y9en3NaCoaf6l9S5
P7/ceWbhqTVifCXlr0I/B0RERAyvFD9IcZE2zpHyQlL8oz1eylMN3/rSqu/N
S7vUkPEjb2mscSXP5Bb1HvGYGF9Jn0s+ChERMfvU5mG0eSTpeOlzteeRfi7d
rxSnSfNxUnwlXb80PyjVX4V+/oiIiPg3pXyIFD9I8ZI2zpGuRzqPFP9I16+9
X+nnf98j78jfvZFxS59cmlcyrdrNnXVZyZNdH1HP90nHS3FX6O8FIiIi/s2F
Wyrbqnu1usq8hj0DitOuev+Y/FUDGtX5E+n82vkyKb7Szq9J16PNj0nzcVZx
mvS5Un176O8LIiJiNivFFdr3snQe7c+t8k7afJH0udp6b22+Szv/KNW3Sz8P
/f1CRETMBqV4QKqr0cY/2rjFKu6yyvNIx5e8sqBwQ1nKTf7BtpOZ47VuUr/X
G4v7lIjjKf1cOr827pLWD4b+fiEiIiZJqU5m+eK6stXdM27OhDXv9Znf5KaU
zTzZtOCAul5IW0etrefRnkc7f6f9uRQvSdepzZtp82DSuGk/FxERsSMo1R3N
+GruNTeVZNybxTu/cNXjLe7ifqmjV+yvNVuPpr1ObT7Hd15LG4dI8ZJVPKaN
V33Ht4iIiNmkFBdJ9S1WfYq08ZVV/bZVXZD2c6X6JW2fKOk5+o6vyDshImJH
UJtH0q4X08Yz2uvXXo9V/CD2jRTqtaR44930+w8fP5Jys3f0+svsinp387TL
Ni05vtX7/KBV3k86XuovGvr7joiI2J7Se02bR5LOb5VvsbovbTzmOw6Rjrda
H6edH7SKM7V9PkP/O0BERGxPqd+RNg+gnbfyPV9jFT9YxQPaOOrHc048+cF1
aTf4pmM5H+w+5H7tXnqh+E/vqvNXoeZDtefRru9DRET0ofSekvJI2v7PVnXR
oc5v9d63ihPuGbhlUsmulLu008rTY/rXu7yq3DkHT1Sa1aVr66+0z8Wqn6dV
vhEREfFclPJF2jptq/e11XvWdxxldby2XnrYpvrZU3ql3e4f/rBTecEh919F
Sx579Ym31fvt+t4v2GqfYqv+D4iIiOeiVKctxUVW67x85zG0WvXxtqqv1s7T
dbnhi6Wj7ki7+V0euODUlMPuwy4r/nPi7P82W+9vFV9px1N7finPabUuABER
O6bSPiNSfskqz+N73k06XttPwPf4a+uLpOOl9XGx9WEI9dylvwt8P19ERMwO
myvyD7w4K+MeundG15xBLa7o+cKcuu616ve19nOt1lv5rp+xOr9VXbc0btJ8
nFX+ymqfF+33QTpemzdjX2BERDwXpfdRz59fWfT01Rn3teLGtTuPNrmFpWcq
vjvioPc8Umxxke+6Kav7kp6jVO8tHW8VD1vNP1r1cZLiKO16wM/77w0REZOp
VPer7Sdp9XOr96l2HJIS/2jjOim/ZHV+q/4J2uO146y9X/oSICJ2TKX3Uaj6
be11hho333GUVf8B6Ty521ZMGDs47TZXPJW7r6DWdZ6x5PD66U94fy5W4+O7
T7v0c2l/58/7fUJExLiV1lNr+2xLat9HVv2Zfc+vaT9Xe/2+97Gt3f5q96F7
U+6pweuq73ux1s2/45Jxt/x4rlmcpn2OvuNt7Thr67joS4CImN1K69qk3/+x
7SciaVU35Xv8Q/W7lvJLVtdj9Vy0z8t3PZW2L6vv7w8iItoq1V34XtcWW32s
7zqr2PokaOfjrPJXvp+jVdyofS7S8VI+lr4EiIjJUsoj+Z7/srr+pM+vSVqt
x9eex/d8nO/9X3z3v7K6fvJOiIhxqq1H8r2+LJSh6pC157fq36j9XCle8r1u
0WrcfM8PauuXpPMTLyEixqm0v4M072b1ub7Xi1nFCb7zUdrrt8oXafMhJa8s
KNxQlnKTf7DtZOZ4rZvU7/XG4j4lZn0mY8vLac9vld9jng4RMazS7/PqgW2j
No7KuC23/Wz0grua3adv7+1cVvOH6NbjS1q9B7V5hqTEA1b76GnjJavrTEp/
dd9xNXknRMTzozTvJtUp+b6eUPNZ0vFWfYR87+vhu47aKl6y6hfqO0/oO76y
el7Sv9NQ+xIiImar0no3aT5OOk/S65QkY6vLsspXWO3Doq1fCtWnNDat8mnS
v1Ptv19ERGxfqS4itjol3/tWSOex+jvd9/vd9/4j0nneTb//8PEjKTd7R6+/
zK6odzdPu2zTkuNbvc/H+e7XFGofHKs4ink6RERbpd+r0jxdbH/X+76e2PZb
0dYdWdVfWc3HJaWuyer5+q6Xk84v/R2UrXlgRETfSnl73+tuYvu9bVWvkpT6
ZEltfkk7H+d7v2Pf/Qd8X6fV90SKn0PVJSIiJl0pvyTNx/l+X4cy1HvZKn/l
u95b+lyrfgJW+9z5/j74zl9Z9aHS5p2sxg0RMVttrsg/8OKsjHvo3hldcwa1
uKLnC3PqutcGW18TKg6xOr/v+nDffSm1x98zcMukkl0pd2mnlafH9K93eVW5
cw6eqPRedxTbvLDv+Uer80j5Jfo4ISK277HXGh49+pOM61F+xUXf2tniypeX
1t1ZXhdsn3TtfJD2/L77PVpdj9V4WuU3pPO0TS0ff/n6lNtaMPRU/5I69+eX
O88sPLXG+/fB6nsSqp+A7+cr/VzKG0u/B3xfPyJiUqTe2/b8sY2P1f2G6u8t
aTUf5/v5+r5fq74QUp75836fEBGzRd/750rG9ve41XsnVH/vUH0+rfpVaj/X
dzwWqs+VdD2+P7fnz68sevrqjPtacePanUeb3MLSMxXfHXEw+O8nRMRYfGDw
xEXN92Vcv5GXdKk+0OwO9djQc83oGu/7n8ZmbHmkUOOsff/6zi+Fqm+3ei6h
6sO1caaUZ6ZPOCLiWZcvritb3T3j5kxY816f+U1uStnMk00LDgR7X1vVLceW
v/L9udL8qe88Ru62FRPGDk67zRVP5e4rqHWdZyw5vH76E8HyPLH1JUhKvCHF
n1L//9DXi4h4vpV+H0r1n77jKOm9n5T3jlbf+9P5rk8evvWlVd+bl3apIeNH
3tJY40qeyS3qPeKxYHXRob4nvvNpVtejPQ/9BxAR23f4oIG/6rsj7fIPNQw8
U9DoWhbvqpv07O+jm4+zWu8W6jpDjZukNs6Z+taStsVVaXfla4On//K5Gtd1
ZUmPx7vOTsw4+N7nN+l5USleiu05IiKG8uoZr674bVvaVQ6b1jz50iY3amHD
mAdrPwq2jj62dVVW9xvb+iyr+bhQ9eehnqN0fGz9prTjJuV7Q/UbQUSMTenv
Smk9XWz9rq2Oj02rdWfa8bGKl2Lra+R7X5LY6vatxod9VRAR29f3PF3S4xmr
+9WOQ6j+D9r6Jav5qVD5Dd915r7rr6ziKKnvpe/xR0RMitI8ne/3dWz1SEm/
Tu31h4qXJEPtM6I9v9X1x1bXxHo6RMT21a6nk86TrXXX2j7MVtcZal8Vbbzk
+3mF6ouV9PjHahyYp0NEbN+krKcLVe8tGdv6vqTkl2J7XlbHW11PqL4K9B9A
RGxfaT8F37+3k/L+DbWPmNV1hsovSVr1n4+tr6bV+Giv0+r8Ul0T6+kQEc8q
/T6U9ltJSpyj1aquJrY8Rmz5pVD7oYTaJy7U89WeR7rOeZtz7p7flHEX7x6+
u6ou5Sp3dRt3benh4P9eERFjseBH87q2LUu7wmHLSvuWNbi5RUMmvvPOnsTH
S77XScV2X7HFS77vN1v3B7RS+32uHtg2auOojNty289GL7ir2X369t7OZTV/
CH4fiIixKOWXtPn5UP15JLXzYrHVe1ut5wrVrzK28bG6/qTXgUtKdeBJ+TsC
EdG3Up2JVO+U9L/rrfo/J6XuJVvzS6Gux+r6Yxt/aT6OuiZExPZdf+uDmz75
OO2OdfvKRRvXNbqez01cWn/9/uDXFbvS+yVUvbe0f1yo/FJscWZshlp/IdWB
068JEbF9Q+23EsqkzAdJx/uej7Man1D9r3wb27yV1T509B9ARPxsJqWPU0fT
Kr8U23OMLU8V2/holZ77wi2VbdW9Wl1lXsOeAcVpV71/TP6qAY1uxldzr7mp
JOPeLN75haseb3EX90sdvWJ/bfD7QESMXWm/Fd/5iqRLvff5yRf57ncaKj8W
qi8TfQYQET+bsc3TdbS4S9J3v8rY9seJ7Xq0x1v187RSuh7yS4iItlIffn6k
3ju79b0eTRpn6XObK/IPvDgr4x66d0bXnEEtruj5wpy67sRLiIifVe1+K9n6
fkx6PwHf92W1fjC2fl/acYvt+y89F+0+3YiI2L5S/YPUNzi290VS7Gj9vWPr
fxXb+Fhdv1SXxXwcIuL50WqeLunvKd/3G9t+u6HqvX3HUdp6cqtx8P39l+JS
qS7R6nMREfGsUn7Jd1/Hjqb0vkt6fkmrNn8S23WGyo9JP5f2Swo1boiIHc1Q
+3smPR6Q1NZ7++4j7dtsfY6+lcZNqlNavriubHX3jJszYc17feY3uSllM082
LTgQ/D4QETuK0t+t2Zr3kLS6r6TUe3c0Y+t7IH5PhD607IeCiBhWqT7cav+F
bH2/+673DnX9occ1dq3qoKR6JKkPbej7RkTE/9+k10vE1k8gVP+lbDUp8Z5U
ryXlkdgnDhExWUq/t6V8VFLeX1qt4pmk7x+XdEOtQ5TqkcgjISJmt9q8U7a+
97V9GqV4yep6snWcfWtVPyaui6QeCRGxQyr1zYttnRfzcckyKfGedr2/1G8/
9H0gImIY6R9+Vqt4KSn7jGSr2v3dCn40r2vbsrQrHLastG9Zg5tbNGTiO+/s
8b4PHSIiJstQfZx8a9Xfm/zS+dF3P3MpX5T09RGIiHh+9L2PQ1LyVEmJl5Iy
nqHGQfreSvsKxTYfjYiIyTK2vJPvOCEp8ZJvY4vHrPpGst4fERF9OG9zzt3z
mzLu4t3Dd1fVpVzlrm7jri097BZuqWyr7tXqKvMa9gwoTrvq/WPyVw2Iry9B
bPvtxjY+vvU9v0bfSEREjFmreo/Y4oek5JdiGzff9yXNl0nza+SREBExBqW/
60PtX+97/7jY4qVQWo2z9jzSvDDr/RERMYlq31+J2UeV+TjTcZaOl9bvS3kk
qV996PtGRERsT2mezncdi5Wx5ZeyNb6Sxk3KI8W27gAREfHzKOUHYnvfafNL
7B/X/rhp71eqL6IPEiIidmSl9d3S8Vb9sWNbH2d1nbGdX1KaL5PiIubXEBGx
IyvtWxdbXwLyS+2rXb8mxUWsX0NERDx3tXknrR1tP5RQcZo03yrFyaHGBxER
MYlq46Wk9/eOLe9k1R8ytvo0RETEbFJb36I9v3ZfsKTkl6zqriWleVIpLmL/
NURExPOvVd7JKp4JlV/yHY9J/UWleFXaH8fqehAREfHcleZ9rOqEtevvFr2x
wG3Kz7jxQ3OW5Y+occ++ecF/VO17NNj4WMVpUn9RafxD3S8iIiKeu77rZLTx
kjZusZq3stpPhH5HiIiI2eeFw6tWX3c64ya3Lfrd6aFpd/Pta298+I4Gsd5G
e36reClU/yhpHKT5Tam/aOjnjIiIiPZa5Z1855ckreq0pbiIdf2IiIgo5Z2k
n0vnCRUvaZXm0Xz3s0JERMTsU1ufE1u8JNW3S3ER+4kgIiKiVim/JMUhoeIl
bb8j1q8hIiKib9ff+uCmTz5Ou2PdvnLRxnWNrudzE5fWX7/fe7wkxW9SXKSd
T0RERES0UsrPSP0YtfGSlL+S5geJixARETEpSv0bpXhJiq+k/BVxESIiIma7
tXPzNoypPfDX/9+pE/2xEREREdv3fwF2N8fK
"], {{0, 0}, {99, 99}}, {0, 1}], Frame -> True, FrameLabel -> {None, None}, FrameTicks -> {{{{98.9806640625, 1}, {60.3474609375, 1000}, {21.6755859375, 2000}, {0.0193359375, 2560}}, {{98.9806640625, 1}, {60.3474609375, 1000}, {21.6755859375, 2000}, {0.0193359375, 2560}}}, {{{0.0193359375, 1}, {38.652539062500004`, 1000}, {77.3244140625, 2000}, {98.9806640625, 2560}}, {{0.0193359375, 1}, {38.652539062500004`, 1000}, {77.3244140625, 2000}, {98.9806640625, 2560}}}}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImageSize -> {615.75, Automatic}, Method -> {"GridLinesInFront" -> True, "AxisPadding" -> Scaled[0.02], "DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultPlotStyle" -> Automatic, "DomainPadding" -> Scaled[0.02], "RangePadding" -> Scaled[0.05]}]Graphics[Raster[CompressedData["
1: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
X30Rc7k=
"], {{0, 0}, {96, 96}}, {0, 1}], Frame -> True, FrameLabel -> {None, None}, FrameTicks -> {{{{95.97500000000001, 1}, {71.025, 500}, {46.025000000000006`, 1000}, {21.025000000000002`, 1500}, {0.025, 1920}}, {{95.97500000000001, 1}, {71.025, 500}, {46.025000000000006`, 1000}, {21.025000000000002`, 1500}, {0.025, 1920}}}, {{{0.025, 1}, {24.975, 500}, {49.975, 1000}, {74.97500000000001, 1500}, {95.97500000000001, 1920}}, {{0.025, 1}, {24.975, 500}, {49.975, 1000}, {74.97500000000001, 1500}, {95.97500000000001, 1920}}}}, GridLinesStyle -> Directive[GrayLevel[0.5, 0.4]], ImageSize -> {597., Automatic}, Method -> {"GridLinesInFront" -> True, "AxisPadding" -> Scaled[0.02], "DefaultBoundaryStyle" -> Automatic, "DefaultGraphicsInteraction" -> {"Version" -> 1.2, "TrackMousePosition" -> {True, False}, "Effects" -> {"Highlight" -> {"ratio" -> 2}, "HighlightPoint" -> {"ratio" -> 2}, "Droplines" -> {"freeformCursorMode" -> True, "placement" -> {"x" -> "All", "y" -> "None"}}}}, "DefaultPlotStyle" -> Automatic, "DomainPadding" -> Scaled[0.02], "RangePadding" -> Scaled[0.05]}]


MatrixPlot[\[CapitalPi]["A"]+ \[CapitalPi]["iA"]]


Aname = "../data/5_4_modulo_2_A.reg";
Bname = "../data/5_4_modulo_2_B.reg";
ABname = "../data/5_4_modulo_2_AB.reg";


Adata = Import[Aname,"Table"];
Bdata = Import[Bname,"Table"];
ABdata = Import[ABname,"Table"];


Aname = "../data/5_4_open_4_A.reg";
Bname = "../data/5_4_open_4_B.reg";
ABname = "../data/5_4_open_4_AB.reg";

Adata = Import[Aname,"Table"];
Bdata = Import[Bname,"Table"];
ABdata = Import[ABname,"Table"];





\[CapitalPi]["A"] = SparseArray[Amat  ,{120,120} ]


s=SparseArray[{{1,1}->1,{2,2}->2,{3,3}->3,{1,3}->4}]







