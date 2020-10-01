(* ::Package:: *)

(* FindConeAngle package for Mathematica, written by Wesley D. Allen*)
(* December 9, 2012 version *)

(* INPUT *)
(* XAtoms = array containing atomic symbols and Cartesian coordinates (in \[CapitalARing]) in the form {{"symbol1",{x1,y1,z1}},{"symbol2",{x2,y2,z2}},...};  If "symbol" is a number rather than an atomic symbol, then that number is used for the corresponding van der Waals radius instead of the default value;
Apex = integer specifying which is the apex atom in XAtoms;
Ligands = integer array in the form {{i1,i2,i3,...},{j1,j2,j3,...},...} specifying which atoms in XAtoms correspond to each ligand for which a cone angle is to be computed;
kPrint = 0  No printing within package;
kPrint = 1  Print {ConeAngle,ConeAxis,ConeAtoms} for each ligand; 
kPrint = 2  For each ligand print both {ConeAngle,ConeAxis,ConeAtoms} and a 3D plot showing the ligand placed inside the cone;
kPrint = 3  For each ligand print both {ConeAngle,ConeAxis,ConeAtoms} and a table of van der Waals radii, vertex angles, and Cartesian coordinates;
kPrint =4  For each ligand print both {ConeAngle,ConeAxis,ConeAtoms} and a histogram of candidate three-atom cone angles (deg) in the range of (\[Theta]2max,\[Theta]cm) *)  

(* OUTPUT *)
(* In addition to any internal printing the package returns an array with {ConeAngle,ConeAxis,ConeAtoms} for each ligand.
ConeAtoms specifies the atoms that are tangent to and fix the cone *)

(*RvdW = Default van der Waals radii (in \[CapitalARing]) from Bondi (1964);
NThreshold = numerical threshold for determining if an atom is contained in a cone;
nLigand = number of atoms in ligand k;
XLigand = array to hold Cartesian coordinates of ligand k;
RvdWLigand = van der Waals radii of atoms in ligand k;
XApex = Cartesian coordinates of apex atom;
rL = distances from apex to ligand atoms;
mL = unit vectors from apex to ligand atoms;
\[Beta]L = vertex angles (rad) of ligand atoms;  *)

FindConeAngle[XAtoms_,Apex_,Ligands_,kPrint_]:=
Module[{ConeAtoms,ConeAxis,ConeAngle,RvdW,nLigand,XLigand,RvdWLigand,XApex,XTangent,kvdW,rL,rCM,nCM,mL,\[Beta]L,\[Beta]12,\[DoubledGamma]L,\[DoubledGamma]pq,\[Alpha],\[Alpha]2max,\[Alpha]cm,\[Alpha]nL,\[Alpha]nL1,\[Alpha]nL2,\[Alpha]sol,\[Alpha]phys,\[Alpha]test,n\[Alpha],naxis,nphys,ntest,tconeL,tpos,NThreshold,Cone1,Cone2,Cone3,Cone3aux,ConeAtomAux,zcone,A123,B123,C123,N123,D123,u123,v123,s123,\[DoubledGamma]123,uN,vN,wN,aw,bw,cw,dw,ax,bx,cx,w,i,j,k,i1,i2,i3,i4,p,q,q1,q2,q3,p1,p2,p3,qC,r,s1,s2,x1,x2,x3,y1,y2,y3,nk1,nk2,nk3,gsphere,Urot,\[Phi],CrossXY,Exclude2Sol,ThreeCone,EndPlot,EndCone,LigandPlot,Ligandtable,AtomColor,ConeResults,TanPts},
RvdW={{"H","He","C","N","O","F","Ne","Si","P","S","Cl","Ar","As","Se","Br","Kr","Te","I","Xe","Fe"},{1.20,1.40,1.70,1.55,1.52,1.47,1.54,2.10,1.80,1.80,1.75,1.88,1.85,1.90,1.85,2.02,2.06,1.98,2.16,2.00},
{White,White,Gray,Blue,Red,Yellow,White,Gray,Orange,Yellow,Green,White,White,White,White,White,White,White,White,Orange}};
Clear[CrossXY,Urot,zcone,gsphere,\[Alpha],\[Phi],x1,x2,x3,y1,y2,y3];
(* Use this explicit function instead of inefficient Mathematica Cross *)
CrossXY[{x1_,x2_,x3_},{y1_,y2_,y3_}]={x2 y3-x3 y2,x3 y1-x1 y3,x1 y2-x2 y1};
(* Unitary matrix for rotation by angle \[Phi] about {x1,x2,x3} axis *)
Urot[\[Phi]_,x1_,x2_,x3_]={{x1^2+(x2^2+x3^2) Cos[\[Phi]],x1 x2 (1-Cos[\[Phi]])+x3 Sin[\[Phi]],x1 x3 (1-Cos[\[Phi]])-x2 Sin[\[Phi]]},{x1 x2 (1-Cos[\[Phi]])-x3 Sin[\[Phi]],x2^2+(x1^2+x3^2) Cos[\[Phi]],x2 x3 (1-Cos[\[Phi]])+x1 Sin[\[Phi]]},{x1 x3 (1-Cos[\[Phi]])+x2 Sin[\[Phi]],x2 x3 (1-Cos[\[Phi]])-x1 Sin[\[Phi]],x3^2+(x1^2+x2^2) Cos[\[Phi]]}};
zcone[\[Alpha]_,x1_,x2_,x3_]=x3^2 (Sin[\[Alpha]])^2-(x1^2+x2^2)(Cos[\[Alpha]])^2; 
gsphere[r_,x1_,y1_,x2_,y2_]=(r^2-(x1-x2)^2-(y1-y2)^2)^(1/2);
NThreshold=10^-8;Exclude2Sol=False;ConeResults={};
XApex=XAtoms[[Apex,2]];
(* Loop over all ligands *)
Do[
nLigand=Length[Ligands[[k]]];
XLigand=Table[XAtoms[[Ligands[[k,j]],2]]-XApex,{j,nLigand}]//N;
(* Collect van der Waals radii of ligand atoms *)
kvdW=Table[Flatten[Position[RvdW,XAtoms[[Ligands[[k,j]],1]]]],{j,nLigand}];
RvdWLigand=Table[XAtoms[[Ligands[[k,j]],1]],{j,nLigand}];
Do[If[Length[kvdW[[j]]]!= 0,If[kvdW[[j,1]]==1,RvdWLigand[[j]]=RvdW[[2,kvdW[[j,2]]]]]],{j,Length[kvdW]}];
(* Vectors and distances from apex to ligand atoms *)
rL=Table[Norm[XLigand[[j]]],{j,nLigand}];
mL=Table[XLigand[[j]]/rL[[j]],{j,nLigand}];
(* Vertex angles (rad) of ligand atoms *)
\[Beta]L=Table[ArcSin[RvdWLigand[[j]]/rL[[j]]],{j,nLigand}];
(* Direction vector to centroid of ligand *)
rCM=Sum[XLigand[[j]],{j,nLigand}]/nLigand;
(* Upper bound for half cone angle *)
If[Abs[Norm[rCM]]>= NThreshold,nCM=rCM/Norm[rCM];\[Alpha]cm=Max[Table[\[Beta]L[[j]]+ArcCos[Dot[mL[[j]],nCM]],{j,nLigand}]],\[Alpha]cm=\[Pi];nCM={0,0,0}];
(* Angles (rad) between direction vectors to ligand atoms *)
\[DoubledGamma]L=Table[0,{i1,nLigand},{i2,nLigand}];
Do[\[DoubledGamma]L[[i1,i2]]=ArcCos[Dot[mL[[i1]],mL[[i2]]]];\[DoubledGamma]L[[i2,i1]]=\[DoubledGamma]L[[i1,i2]],{i1,2,nLigand},{i2,i1-1}];
If[(kPrint== 3),Print["Ligand atom van der Waals radii (r), vertex angles (\[Beta], deg), and Cartesian coordinates (x,y,z)"];
Ligandtable=Table[Flatten[{XAtoms[[Ligands[[k,j]],1]],RvdWLigand[[j]],\[Beta]L[[j]]/Degree,XAtoms[[Ligands[[k,j]],2]]}],{j,nLigand}];
Print[TableForm[Join[{{" ","  r","  \[Beta]","  x","  y","  z"}},Ligandtable]]]];
ConeAtoms={};XTangent={};
(* tconeL(i,j) test of whether atom i lies completely within cone of atom j *)
tconeL=Table[\[Beta]L[[i2]]>=\[Beta]L[[i1]]+\[DoubledGamma]L[[i1,i2]] ,{i1,nLigand},{i2,nLigand}];
Do[tconeL[[i1,i1]]=False,{i1,nLigand}];
tpos=Position[tconeL,True];
(* Ligand atoms that lie in shadow of another ligand atom *)
qC=Union[Table[tpos[[i1,1]],{i1,Length[tpos]}]];
(* Active single candidates for minimum cone *)
Cone1=Complement[Range[nLigand],qC];
nk1=Length[Cone1];
If[nk1== 1,ConeAtoms=Cone1;
j=Cone1[[1]];ConeAxis=mL[[j]];\[Alpha]=\[Beta]L[[j]];ConeAngle=2 \[Alpha]/Degree;ConeAtoms[[1]]=Ligands[[k,j]];Goto[EndCone]];
If[Exclude2Sol,Print["Two atom solutions excluded"];\[Alpha]=Max[\[Beta]L];\[Alpha]cm=\[Pi];Goto[ThreeCone]];
(* Active candidates pairs for minimum cone *)
Cone2=Flatten[Table[{Cone1[[i1]],Cone1[[i2]]},{i1,1,nk1},{i2,i1+1,nk1}],1];
nk2=Length[Cone2];
(* Cones angles for active candidates pairs *)
\[Alpha]nL=Table[{i1,i2}=Cone2[[p]];(\[Beta]L[[i1]]+\[Beta]L[[i2]]+\[DoubledGamma]L[[i1,i2]])/2,{p,nk2}];
(* Select max cone angle (\[Alpha]) and compute cone axis for the \[Alpha] pair *)
\[Alpha]=Max[\[Alpha]nL];{i1,i2}=Cone2[[Position[\[Alpha]nL,\[Alpha]][[1,1]]]];\[Beta]12=\[DoubledGamma]L[[i1,i2]];
(* Compute cone axis test if all atoms are contained completely inside cone for the \[Alpha] pair *)
If[Chop[Sin[\[Beta]12],NThreshold]>0,naxis=Csc[\[Beta]12](Sin[(\[Beta]12+\[Beta]L[[i1]]-\[Beta]L[[i2]])/2] mL[[i1]]+Sin[(\[Beta]12-\[Beta]L[[i1]]+\[Beta]L[[i2]])/2]mL[[i2]]);
naxis=naxis/Norm[naxis];
\[DoubledGamma]pq=Table[ArcCos[Dot[mL[[Cone1[[p]]]],naxis]],{p,nk1}];
tconeL=Table[\[Alpha]>= \[DoubledGamma]pq[[p]]+\[Beta]L[[Cone1[[p]]]],{p,nk1}],\[Alpha]=Max[\[Beta]L];\[Alpha]cm=\[Pi];Goto[ThreeCone]];
If[Count[tconeL,True]== nk1,ConeAngle=2 \[Alpha]/Degree; ConeAxis=naxis;
ConeAtoms={Ligands[[k,i1]],Ligands[[k,i2]]};
XTangent=Table[i=Part[{i1,i2},p];
RvdWLigand[[i]]Cot[\[Beta]L[[i]]]Csc[\[Alpha]-\[Beta]L[[i]]](Sin[\[Alpha]]mL[[i]]-Sin[\[Beta]L[[i]]] naxis),{p,2}];Goto[EndCone]];
Label[ThreeCone];
(* Construct candidate triples *)
Cone3=Flatten[Table[{Cone1[[i1]],Cone1[[i2]],Cone1[[i3]]},{i1,nk1},{i2,i1+1,nk1},{i3,i2+1,nk1}],2];nk3=Length[Cone3];
\[Alpha]2max=\[Alpha]  (* Save best cone angle from pair analysis as a lower bound. *);
(* Compute cone angles and axes of all active candidate triples *)
\[Alpha]nL=Table[{i1,i2,i3}=Cone3[[p]];
u123={Cos[\[Beta]L[[i1]]],Cos[\[Beta]L[[i2]]],Cos[\[Beta]L[[i3]]]};
v123={Sin[\[Beta]L[[i1]]],Sin[\[Beta]L[[i2]]],Sin[\[Beta]L[[i3]]]};
N123=Transpose[{CrossXY[mL[[i2]],mL[[i3]]],CrossXY[mL[[i3]],mL[[i1]]],CrossXY[mL[[i1]],mL[[i2]]]}];
uN=Dot[N123,u123];vN=Dot[N123,v123];
\[DoubledGamma]123=Dot[mL[[i1]],CrossXY[mL[[i2]],mL[[i3]]]];
If[Abs[\[DoubledGamma]123]>NThreshold,s123=Sign[\[DoubledGamma]123]{1,1,1,1},s123={1,-1,1,-1}];
A123=Dot[uN,uN];B123=Dot[vN,vN];C123=Dot[uN,vN];D123=\[DoubledGamma]123^2;
ax=A123-B123;bx=A123+B123-2 D123;cx=(2 C123)^2;
aw=cx+ax^2;bw=(ax bx)/aw;cw=(bx^2-cx)/aw;dw=Sqrt[bw^2-cw];
If[Abs[dw]<10 NThreshold,dw=0];
{ax,bx}={ArcCos[-bw-dw]/2,ArcCos[-bw+dw]/2};
\[Alpha]sol={ax,bx,\[Pi]-ax,\[Pi]-bx};
ntest=Table[cx=Cos[\[Alpha]sol[[j]]];bx=Sin[\[Alpha]sol[[j]]];Chop[A123 cx^2+B123 bx^2+2 C123 cx bx-D123,NThreshold]== 0,{j,4}];
qC=Flatten[Position[ntest,True]];
\[Alpha]phys=Table[\[Alpha]sol[[qC[[j]]]],{j,Length[qC]}];
\[Alpha]test=Table[(\[Alpha]phys[[j]]>= \[Alpha]2max)&&(\[Alpha]phys[[j]]<= \[Alpha]cm),{j,Length[\[Alpha]phys]}];
qC=Flatten[Position[\[Alpha]test,True]];
\[Alpha]phys=Sort[Table[\[Alpha]phys[[qC[[j]]]],{j,Length[qC]}],Greater];
n\[Alpha]=Length[\[Alpha]phys];nphys={};
If[n\[Alpha]!= 0,
ax={Abs[Cos[\[DoubledGamma]L[[i1,i2]]]],Abs[Cos[\[DoubledGamma]L[[i2,i3]]]],Abs[Cos[\[DoubledGamma]L[[i3,i1]]]]};
bx={{i1,i2,i3},{i2,i3,i1},{i3,i1,i2}};
{p1,p2,p3}=bx[[Position[ax,Min[ax]][[1,1]]]];
nphys=Table[wN=(Cos[\[Alpha]phys[[j]]-\[Beta]L[[p1]]]-Cos[\[Alpha]phys[[j]]-\[Beta]L[[p2]]] Cos[\[DoubledGamma]L[[p1,p2]]])mL[[p1]]+
(Cos[\[Alpha]phys[[j]]-\[Beta]L[[p2]]]-Cos[\[Alpha]phys[[j]]-\[Beta]L[[p1]]] Cos[\[DoubledGamma]L[[p1,p2]]])mL[[p2]];
cx=Sqrt[1-(Csc[\[DoubledGamma]L[[p1,p2]]])^4 Dot[wN,wN]];
dw=CrossXY[mL[[p1]],mL[[p2]]];
aw=s123[[j]]Sign[Dot[(uN Cos[\[Alpha]phys[[j]]]+vN Sin[\[Alpha]phys[[j]]]),dw]];
(Csc[\[DoubledGamma]L[[p1,p2]]])^2 (wN+cx aw Sin[\[DoubledGamma]L[[p1,p2]]] dw),{j,n\[Alpha]}]];{\[Alpha]phys,nphys},{p,nk3}];
(* Eliminate cases whose cone angles did not fall within lower and upper bounds *)
qC={};Do[If[Length[\[Alpha]nL[[q,1]]]!=  0,qC=Join[qC,{q}]],{q,nk3}];
\[Alpha]nL=Table[\[Alpha]nL[[qC[[q]]]],{q,Length[qC]}];
Cone3=Table[Cone3[[qC[[q]]]],{q,Length[qC]}];
nk3=Length[Cone3];
\[Alpha]nL1=Table[{\[Alpha]nL[[q,1,1]],\[Alpha]nL[[q,2,1]]},{q,nk3}];
(* Add cases with multiple, valid cone angle solutions *)
qC={};Do[If[Length[\[Alpha]nL[[q,1]]]== 2,qC=Join[qC,{q}]],{q,nk3}];
(* Save triples with auxiliary solutions. *)
ConeAtomAux=Table[{i1,i2,i3}=Cone3[[qC[[q]]]];{2 \[Alpha]nL[[qC[[q]],1]]/Degree,{Ligands[[k,i1]],Ligands[[k,i2]],Ligands[[k,i3]]}},{q,Length[qC]}];
(* Join all three-atom solutions *)
\[Alpha]nL2=Table[{\[Alpha]nL[[qC[[q]],1,2]],\[Alpha]nL[[qC[[q]],2,2]]},{q,Length[qC]}];
\[Alpha]nL=Join[\[Alpha]nL1,\[Alpha]nL2];
(* Form histogram of candidate three-atom cone angles *)
If[kPrint== 4,Print["Allowed range of candidate three-atom cone angles = ",2 {\[Alpha]2max,\[Alpha]cm}/Degree];
Print[Histogram[2 Transpose[\[Alpha]nL][[1]]/Degree,10,Frame-> True]]];
Cone3aux=Table[Cone3[[qC[[q]]]],{q,Length[qC]}];
Cone3=Join[Cone3,Cone3aux];nk3=Length[Cone3];
If[kPrint>4,Print["Ligand atoms with auxiliary solutions = ",ConeAtomAux]];
(* tconeL(q,p) tests whether atom Cone1[[p]] lies completely inside the cone for triple Cone3[[q]] *)
\[DoubledGamma]pq=Table[ArcCos[Dot[mL[[Cone1[[p]]]],\[Alpha]nL[[q,2]]]],{p,nk1},{q,nk3}];
tconeL=Table[Chop[\[Alpha]nL[[q,1]]- \[DoubledGamma]pq[[p,q]]-\[Beta]L[[Cone1[[p]]]],NThreshold]>= 0,{q,nk3},{p,nk1}];
(* Find all remaining cones of triples that contain all ligand atoms *)
qC=Flatten[Table[If[Count[tconeL[[q]],True]== nk1,q,{}],{q,nk3}]];
(* Determine minimum cone angle *)
If[Length[qC]>0,\[Alpha]=Min[Table[\[Alpha]nL[[qC[[q]],1]],{q,Length[qC]}]];j=Position[\[Alpha]nL,\[Alpha]][[1,1]];naxis=\[Alpha]nL[[j,2]];{i1,i2,i3}=Cone3[[j]];
XTangent=Table[i=Cone3[[j,p]];
RvdWLigand[[i]]Cot[\[Beta]L[[i]]]Csc[\[Alpha]-\[Beta]L[[i]]](Sin[\[Alpha]]mL[[i]]-Sin[\[Beta]L[[i]]] naxis),{p,3}],\[Alpha]=\[Alpha]cm;naxis=nCM;XTangent={};Print["No tangent solutions found."]];
ConeAngle=2 \[Alpha]/Degree;ConeAxis=naxis;ConeAtoms={Ligands[[k,i1]],Ligands[[k,i2]],Ligands[[k,i3]]};
Label[EndCone];
If[kPrint> 0,Print["Ligand = ",k,"   Ligand atoms forming cone = ",ConeAtoms];
Print["Cone angle (deg) = ",ConeAngle,"   Cone axis = ",ConeAxis]];
(* Plot ligand inside cone *)
If[(kPrint!= 2),Goto[EndPlot]];
(* Rotate ligand atoms to make the cone axis the z-axis *)
{x1,x2,x3}=ConeAxis;
If[Chop[Norm[{x1,x2,x3}]-1,NThreshold]!=  0,Goto[EndPlot]];
If[Chop[Norm[{x1,x2}],NThreshold]!=0,{y1,y2}={x2,-x1}/Sqrt[x1^2+x2^2];XLigand=Table[Dot[Urot@@{-ArcCos[x3],y1,y2,0},XLigand[[j]]],{j,nLigand}];
XTangent=Table[Dot[Urot@@{-ArcCos[x3],y1,y2,0},XTangent[[j]]],{j,Length[XTangent]}]];
AtomColor=Table[Transparent,{j,nLigand}];
Do[If[Length[kvdW[[j]]]!= 0,If[kvdW[[j,1]]==1,AtomColor[[j]]=RvdW[[3,kvdW[[j,2]]]]]],{j,nLigand}];
Do[r=RvdWLigand[[j]];{y1,y2,y3}=XLigand[[j]];
s1[j]=Plot3D[y3+gsphere[r,x1,x2,y1,y2],{x1,y1-r,y1+r},{x2,y2-r,y2+r},PlotStyle-> AtomColor[[j]],Mesh->None,Lighting->"Neutral"];
s2[j]=Plot3D[y3-gsphere[r,x1,x2,y1,y2],{x1,y1-r,y1+r},{x2,y2-r,y2+r},PlotStyle-> AtomColor[[j]],Mesh->None,Lighting->"Neutral"],{j,nLigand}];
If[Length[ConeAtoms]>1,TanPts={Graphics3D[Join[{Black},Table[Sphere[XTangent[[j]],0.15],{j,Length[XTangent]}]]]},
TanPts={ParametricPlot3D[RvdWLigand[[Cone1[[1]]]]Cos[\[Alpha]]{Cos[\[Phi]],Sin[\[Phi]],Cot[\[Alpha]]},{\[Phi],0,2\[Pi]},PlotStyle-> Thick]}];
w=(6/5)Max[Table[rL[[j]]+RvdWLigand[[j]],{j,nLigand}]];
LigandPlot=Show@@Join[{ContourPlot3D[zcone[\[Alpha],x1,x2,x3]== 0,{x1,-w,w},{x2,-w,w},{x3,0,w Sign[Cos[\[Alpha]]]},PlotPoints-> 30,ColorFunction->"Aquamarine",Lighting->"Neutral"]},TanPts,Flatten[Table[{s1[j],s2[j]},{j,Length[XLigand]}]],{ImageSize-> 500,BoxRatios->{1,1,1},Boxed-> False,Axes-> False,PlotRange->{{-w,w},{-w,w},{-w,w}}}];
Print[LigandPlot];
Label[EndPlot];
ConeResults=Join[ConeResults,{{ConeAngle,ConeAxis,ConeAtoms}}],{k,Length[Ligands]}];ConeResults];
