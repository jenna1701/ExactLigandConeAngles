(* ::Package:: *)

ligandsIn=OpenRead["ComplexDataBase1.txt"];
LigandTag=StringJoin["(X",ToString[ComplexStructure],")"];
readtest=True;
xC={};
While[readtest,xr=Read[ligandsIn,Word];
If[xr== EndOfFile,readtest=False];
If[xr== LigandTag,
ligandname =Read[ligandsIn,Word];
While[readtest,k=Read[ligandsIn,Number];
If[k== 0,Break[]];rC=Read[ligandsIn,{Word,Number,Number,Number}];
xC=Join[xC,{{rC[[1]],{rC[[2]],rC[[3]],rC[[4]]}}}]]]];
Apex=1;
Ligands={Range[2,Length[xC]]};
XAtoms=xC;
Print[ligandname ];
