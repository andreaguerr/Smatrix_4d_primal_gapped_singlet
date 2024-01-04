(* ::Package:: *)

(* ::Subsection::Closed:: *)
(*Basic constants and parameters*)


(*Set your home directory*)
ClearAll[dir];
dir=SetDirectory[NotebookDirectory[]];


(* precision parameters set to a generous value *)
Unprotect[precgrid,prec]; precgrid= 350; prec=250;Protect[precgrid,prec];
$MaxExtraPrecision=100000;


(* Loading SDPB auxiliary file*)
<< (dir<>"/SDPB_memoized.m")


(* Set your favorite directory to store the integrals of partial waves *)
ClearAll[integralsdir];
integralsdir=dir<>"/integral_storage_provv";


(* Choice of the Ansatz parameters, minimal choice set. *)
ClearAll[\[Sigma]center0,\[Sigma]centerlist,lengthcenters];
\[Sigma]center0=20/3;
\[Sigma]centerlist={};
lengthcenters=0;


(* ::Subsection::Closed:: *)
(*Variables and Grids*)


(* Mandelstam and conformal mappings *)
ClearAll[t,u,\[Rho]\[Rho],\[ScriptS]];
t[s_,x_]:=1/2 (-4+s) (-1+x); 
u[s_,x_]:=-(1/2) (-4+s) (1+x);
\[Rho]\[Rho][s_,\[Sigma]_]:=(Sqrt[\[Sigma]-4]-Sqrt[4-s])/(Sqrt[\[Sigma]-4]+Sqrt[4-s]);
\[ScriptS][\[Rho]_,\[Sigma]_]:=-((-8+\[Rho]^2 (-8+\[Sigma])+\[Sigma]-2 \[Rho] \[Sigma])/(1+\[Rho])^2);


(* Choice of the grids *)
ClearAll[npts0,npts1,npts,\[Phi]gridCheb0,\[Phi]gridCheb1,\[Rho]gridCheb0,\[Rho]gridCheb1,sgrid];
npts0=250;npts1=50;
npts=npts0+npts1*lengthcenters;
\[Phi]gridCheb0=SetPrecision[\[Pi]/2 (1+Cos[\[Pi] Range[npts0+1]/(npts0+1)])//Drop[#,-1]&,precgrid]//Reverse;
\[Phi]gridCheb1=SetPrecision[\[Pi]/2 (1+Cos[\[Pi] Range[npts1+1]/(npts1+1)])//Drop[#,-1]&,precgrid]//Reverse;
\[Rho]gridCheb0=SetPrecision[(Exp[I \[Phi]gridCheb0]),precgrid];
\[Rho]gridCheb1=SetPrecision[(Exp[I \[Phi]gridCheb1]),precgrid];
sgrid=Chop[\[ScriptS][\[Rho]gridCheb0,\[Sigma]center0],10^-prec]~Join~Flatten[Table[Chop[\[ScriptS][\[Rho]gridCheb1,\[Sigma]],10^-prec],{\[Sigma],\[Sigma]centerlist}],1];


(* ::Subsection::Closed:: *)
(*Amplitude Ansatz*)


(* Crossing symmetrized base monomials *)
ClearAll[\[Rho]disc,\[Rho]ddisc];
\[Rho]disc[s_,t_,u_][\[Sigma]_,n_]:=\[Rho]\[Rho][s,\[Sigma]]^n+\[Rho]\[Rho][t,\[Sigma]]^n+\[Rho]\[Rho][u,\[Sigma]]^n;
\[Rho]ddisc[s_,t_,u_][\[Sigma]_,n_,m_]:=\[Rho]\[Rho][s,\[Sigma]]^n (\[Rho]\[Rho][t,\[Sigma]]^m+\[Rho]\[Rho][u,\[Sigma]]^m)+\[Rho]\[Rho][s,\[Sigma]]^m (\[Rho]\[Rho][t,\[Sigma]]^n+\[Rho]\[Rho][u,\[Sigma]]^n)+(\[Rho]\[Rho][u,\[Sigma]]^n \[Rho]\[Rho][t,\[Sigma]]^m+\[Rho]\[Rho][u,\[Sigma]]^m \[Rho]\[Rho][t,\[Sigma]]^n);


(* Free variables *)
ClearAll[\[Alpha]abc,\[Alpha],\[Beta]];
\[Alpha]abc[Nmax_,Mmax_]:={th}~Join~Table[\[Alpha][0][n],{n,0,Nmax}]~Join~Flatten[Table[Table[\[Alpha][\[ScriptI]][n],{n,1,Mmax}],{\[ScriptI],1,lengthcenters}]]~Join~Flatten[Table[If[n+m<=Nmax,\[Beta][0][n,m],Nothing],{n,1,Nmax},{m,1,n}]]~Join~Flatten[Table[Table[If[n+m<=Mmax,\[Beta][\[ScriptI]][n,m],Nothing],{n,1,Mmax},{m,1,n}],{\[ScriptI],1,lengthcenters}]];


(* Ansatz *)
ClearAll[Amplitude];
Amplitude[s_,t_,u_][Nmax_,Mmax_]:={1/(\[Rho]\[Rho][s,\[Sigma]center0]-1)+1/(\[Rho]\[Rho][t,\[Sigma]center0]-1)+1/(\[Rho]\[Rho][u,\[Sigma]center0]-1)}~Join~Join[Table[\[Rho]disc[s,t,u][\[Sigma]center0,n],{n,0,Nmax}],Flatten[Table[\[Rho]disc[s,t,u][\[Sigma],n],{n,1,Mmax},{\[Sigma],\[Sigma]centerlist}]],Flatten[Table[If[n+m<=Nmax,\[Rho]ddisc[s,t,u][\[Sigma]center0,n,m],Nothing],{n,1,Nmax},{m,1,n}]],Flatten[Table[If[n+m<=Mmax,\[Rho]ddisc[s,t,u][\[Sigma],n,m],Nothing],{n,1,Mmax},{m,1,n},{\[Sigma],\[Sigma]centerlist}]]];


(* ::Subsection::Closed:: *)
(*Partial Wave Integrals*)


(* ::Subsubsection::Closed:: *)
(*Numerical Integrals Generation Section*)


Unprotect[\[ScriptP]recgrid]; \[ScriptP]recgrid=1000;Protect[\[ScriptP]recgrid];
ClearAll[\[ScriptP]goal];\[ScriptP]goal=250;


ClearAll[PWN];
PWN[l_,\[Sigma]_,n_][\[DoubleStruckS]_]:=Block[{prec0,prec1,iter0,iter1,x},
prec0=\[ScriptP]goal;
iter0= (\[Rho]\[Rho][t[\[DoubleStruckS],x],\[Sigma]]^n)LegendreP[l,x]//NIntegrate[#,{x,-1,1},WorkingPrecision->prec0,MaxRecursion->150,AccuracyGoal->prec0/2,PrecisionGoal->prec0/2,Method->{"GlobalAdaptive",Method->"ClenshawCurtisRule"}]&//Quiet;
prec1=prec0+40;
iter1= (\[Rho]\[Rho][t[\[DoubleStruckS],x],\[Sigma]]^n)LegendreP[l,x]//NIntegrate[#,{x,-1,1},WorkingPrecision->prec1,MaxRecursion->150,AccuracyGoal->prec1/2,PrecisionGoal->prec1/2,Method->{"GlobalAdaptive",Method->"ClenshawCurtisRule"}]&//Quiet;
While[Abs[iter0-iter1]>10^-\[ScriptP]goal||\[ScriptP]goal>1500,
prec0=prec1;
prec1=prec0+40;
iter0=iter1;
iter1= (\[Rho]\[Rho][t[\[DoubleStruckS],x],\[Sigma]]^n)LegendreP[l,x]//NIntegrate[#,{x,-1,1},WorkingPrecision->prec1,MaxRecursion->150,AccuracyGoal->prec1/2,PrecisionGoal->prec1/2,Method->{"GlobalAdaptive",Method->"ClenshawCurtisRule"}]&//Quiet;
Print[prec1];];iter1
];

PWN[l_,\[Sigma]_,n_,m_][\[DoubleStruckS]_]:=Block[{prec0,prec1,iter0,iter1,x},
prec0=\[ScriptP]goal;
iter0= (\[Rho]\[Rho][t[\[DoubleStruckS],x],\[Sigma]]^n \[Rho]\[Rho][u[\[DoubleStruckS],x],\[Sigma]]^m)LegendreP[l,x]//NIntegrate[#,{x,-1,1},WorkingPrecision->prec0,MaxRecursion->150,AccuracyGoal->prec0/2,PrecisionGoal->prec0/2,Method->{"GlobalAdaptive",Method->"ClenshawCurtisRule"}]&//Quiet;
prec1=prec0+40;
iter1=(\[Rho]\[Rho][t[\[DoubleStruckS],x],\[Sigma]]^n \[Rho]\[Rho][u[\[DoubleStruckS],x],\[Sigma]]^m)LegendreP[l,x]//NIntegrate[#,{x,-1,1},WorkingPrecision->prec1,MaxRecursion->150,AccuracyGoal->prec1/2,PrecisionGoal->prec1/2,Method->{"GlobalAdaptive",Method->"ClenshawCurtisRule"}]&//Quiet;
While[Abs[iter0-iter1]>10^-\[ScriptP]goal||\[ScriptP]goal>1500,
prec0=prec1;
prec1=prec0+40;
iter0=iter1;
iter1=(\[Rho]\[Rho][t[\[DoubleStruckS],x],\[Sigma]]^n \[Rho]\[Rho][u[\[DoubleStruckS],x],\[Sigma]]^m)LegendreP[l,x]//NIntegrate[#,{x,-1,1},WorkingPrecision->prec1,MaxRecursion->150,AccuracyGoal->prec1/2,PrecisionGoal->prec1/2,Method->{"GlobalAdaptive",Method->"ClenshawCurtisRule"}]&//Quiet;
Print[prec1];];iter1
];


ClearAll[PWNList];
PWNList[l_][\[Sigma]_,n_]:=PWNList[l][\[Sigma],n]=ParallelTable[(*Print[ss//N];*)PWN[l,\[Sigma],n][ss],{ss,SetPrecision[sgrid,\[ScriptP]recgrid]}];
PWNList[l_][\[Sigma]_,n_,m_]:=PWNList[l][\[Sigma],n,m]=ParallelTable[(*Print[ss//N];*)PWN[l,\[Sigma],n,m][ss],{ss,SetPrecision[sgrid,\[ScriptP]recgrid]}];


ClearAll[PWNSave];
PWNSave[l_][\[Sigma]_,n_]:=PWNSave[l][\[Sigma],n]=Block[{namefile},
namefile=If[lengthcenters>0,
"rhorhoint_4d_gapped_"<>StringDrop[StringJoin@@Table[ToString[N[\[Sigma]\[Sigma],4]]<>"_",{\[Sigma]\[Sigma],\[Sigma]centerlist}],-1]<>"_"<>ToString[l]<>"_"<>ToString[N[\[Sigma],4]]<>"_"<>ToString[n]<>".txt","rhorhoint_4d_gapped_"<>ToString[l]<>"_"<>ToString[N[\[Sigma],4]]<>"_"<>ToString[n]<>".txt"];
If[FileExistsQ[integralsdir<>"/"<>namefile],Print["Chirp!_"<>namefile],Export[integralsdir<>"/"<>namefile,ToString[FullForm[PWNList[l][\[Sigma],n]]]]]];

PWNSave[l_][\[Sigma]_,n_,m_]:=PWNSave[l][\[Sigma],n,m]=Block[{namefile},
namefile=If[lengthcenters>0,
"rhorhoint_4d_gapped_"<>StringDrop[StringJoin@@Table[ToString[N[\[Sigma]\[Sigma],4]]<>"_",{\[Sigma]\[Sigma],\[Sigma]centerlist}],-1]<>"_"<>ToString[l]<>"_"<>ToString[N[\[Sigma],4]]<>"_"<>ToString[n]<>"_"<>ToString[m]<>".txt","rhorhoint_4d_gapped_"<>ToString[l]<>"_"<>ToString[N[\[Sigma],4]]<>"_"<>ToString[n]<>"_"<>ToString[m]<>".txt"];
If[FileExistsQ[integralsdir<>"/"<>namefile],Print["Chirp!_"<>namefile],Export[integralsdir<>"/"<>namefile,ToString[FullForm[PWNList[l][\[Sigma],n,m]]]]]];


(* ::Subsubsection::Closed:: *)
(*Construction partial wave projection of the ansatz*)


(* Loading seed integrals *)
ClearAll[preintReg];
preintReg[l_][\[Sigma]_,n_]:=preintReg[\[Sigma],n]=
Block[{namefile},
namefile=If[lengthcenters>0,
"rhorhoint_4d_gapped_"<>StringDrop[StringJoin@@Table[ToString[N[\[Sigma]\[Sigma],4]]<>"_",{\[Sigma]\[Sigma],\[Sigma]centerlist}],-1]<>"_"<>ToString[l]<>"_"<>ToString[N[\[Sigma],4]]<>"_"<>ToString[n]<>".txt","rhorhoint_4d_gapped_"<>ToString[l]<>"_"<>ToString[N[\[Sigma],4]]<>"_"<>ToString[n]<>".txt"];
If[FileExistsQ[integralsdir<>"/"<>namefile],ToExpression[Import[integralsdir<>"/"<>namefile]],PWNSave[l][\[Sigma],n];ToExpression[Import[integralsdir<>"/"<>namefile]]]];
preintReg[l_][\[Sigma]_,n_,m_]:=preintReg[l][\[Sigma],n,m]=
Block[{namefile},
namefile=If[lengthcenters>0,
"rhorhoint_4d_gapped_"<>StringDrop[StringJoin@@Table[ToString[N[\[Sigma]\[Sigma],4]]<>"_",{\[Sigma]\[Sigma],\[Sigma]centerlist}],-1]<>"_"<>ToString[l]<>"_"<>ToString[N[\[Sigma],4]]<>"_"<>ToString[n]<>"_"<>ToString[m]<>".txt","rhorhoint_4d_gapped_"<>ToString[l]<>"_"<>ToString[N[\[Sigma],4]]<>"_"<>ToString[n]<>"_"<>ToString[m]<>".txt"];
If[FileExistsQ[integralsdir<>"/"<>namefile],ToExpression[Import[integralsdir<>"/"<>namefile]],PWNSave[l][\[Sigma],n,m];ToExpression[Import[integralsdir<>"/"<>namefile]]]];


(* Building partial wave projection of the base monomials *)
ClearAll[TpwTh,intReg];
intReg[l_][\[Sigma]_,n_]:=1/(32 \[Pi]) (Conjugate[\[Rho]\[Rho][sgrid,\[Sigma]]]^n 2 KroneckerDelta[l,0]+preintReg[l][\[Sigma],n]+(-1)^l preintReg[l][\[Sigma],n]);
intReg[l_][\[Sigma]_,n_,m_]:=1/(32 \[Pi]) (Conjugate[\[Rho]\[Rho][sgrid,\[Sigma]]]^n (preintReg[l][\[Sigma],m]+(-1)^l preintReg[l][\[Sigma],m])+Conjugate[\[Rho]\[Rho][sgrid,\[Sigma]]]^m (preintReg[l][\[Sigma],n]+(-1)^l preintReg[l][\[Sigma],n])+(preintReg[l][\[Sigma],n,m]+(-1)^l preintReg[l][\[Sigma],n,m]));
TpwTh[\[ScriptL]_]:=1/(32 \[Pi]) (-2 KroneckerDelta[0,\[ScriptL]]+(2 KroneckerDelta[0,\[ScriptL]])/(-1+Conjugate[\[Rho]\[Rho][sgrid,\[Sigma]center0]])-(8 Sqrt[2/3] (-2+Sqrt[sgrid])^\[ScriptL] (2+Sqrt[sgrid])^(-1-\[ScriptL]))/(1+2 \[ScriptL]));


(* Building partial wave projection of Ansatz *)
ClearAll[Tpw];
Tpw[\[ScriptL]_][Nmax_,Mmax_]:=Tpw[\[ScriptL]][Nmax,Mmax]=Join[{TpwTh[\[ScriptL]]},Table[intReg[\[ScriptL]][\[Sigma]center0,n],{n,0,Nmax}],Flatten[Table[intReg[\[ScriptL]][\[Sigma],n],{n,1,Mmax},{\[Sigma],\[Sigma]centerlist}],1],Flatten[Table[If[n+m<=Nmax,intReg[\[ScriptL]][\[Sigma]center0,n,m],Nothing],{n,1,Nmax},{m,1,n}],1],Flatten[Table[If[n+m<=Mmax,intReg[\[ScriptL]][\[Sigma],n,m],Nothing],{n,1,Mmax},{m,1,n},{\[Sigma],\[Sigma]centerlist}],1]]//Transpose//Chop[#,10^-prec]&;


(* ::Subsection::Closed:: *)
(*"Improved" positivity constraints for 0<t<4]*)


(* Building improved positivity constraints *)
ClearAll[PositivityList];
PositivityList[Nd_,Md_][tt_]:=PositivityConstraintsSub[Nd,Md][tt]=Block[{t},
Table[Amplitude[s,t,4-s-t][Nd,Md]/.t->tt//Conjugate//Im,{s,sgrid}]//Chop[#,10^-prec]&];
ClearAll[\[ScriptCapitalT]list];
\[ScriptCapitalT]list=Reverse[SetPrecision[2 (1+Cos[\[Pi] Range[10+1]/(10+1)]),precgrid]]~Join~Cases[Reverse[SetPrecision[2 (1+Cos[\[Pi] Range[2000+1]/(2000+1)]),precgrid]][[1;;-1;;2]],a_/;a>3+999/1000]//Sort;
ClearAll[ImprovedPositivityList];
ImprovedPositivityList[Nd_,Md_,Lmax_][tt_]:=ImprovedPositivityList[Nd,Md,Lmax][tt]=Block[{t},PositivityList[Nd,Md][tt]-Sum[16 \[Pi] (2\[ScriptL]+1)LegendreP[\[ScriptL],1+2 t/(sgrid-4)]Im[Tpw[\[ScriptL]][Nd,Md]]/.t->tt,{\[ScriptL],0,Lmax,2}]//Chop[#,10^-prec]&];


(* ::Subsection::Closed:: *)
(*Objective Definition*)


(* Defining the quartic coupling objective *)
ClearAll[\[DoubleStruckCapitalC]0];
\[DoubleStruckCapitalC]0[Nmax_,Mmax_]:=Block[{s,t},1/(32 \[Pi]) Amplitude[s,t,4-s-t][Nmax,Mmax]/.s->4/3/.t->4/3];


(* ::Subsection::Closed:: *)
(*Setup of the Quadratic Program*)


(* Build the quadratic program *)
ClearAll[almostzeroarray];
almostzeroarray[len_,repl_] :=ReplacePart[ConstantArray[0,len],repl];
ClearAll[buildSmatrixProgram];
buildSmatrixProgram[fname_,Nd_,Md_,Lmax_]:=Block[{obj,norm,\[Alpha]vecs,\[Beta]vecs,quadraticconstraints,linearconstraints,tomatrix,idmatrix},

tomatrix[a_]:={{-(1/2)Im[a],Re[a]},{Re[a],2 Im[a]}};
idmatrix:={{1,0},{0,0}};

\[Alpha]vecs=Table[Sqrt[(sgrid-4)/sgrid]Tpw[\[ScriptL]][Nd,Md],{\[ScriptL],0,Lmax,2}]//Flatten[#,1]&;
\[Beta]vecs=ParallelTable[Print[tt//N];ImprovedPositivityList[Nd,Md,Lmax][tt],{tt,\[ScriptCapitalT]list}]//Flatten[#,1]&;

quadraticconstraints =
Join[Table[
PositiveMatrixWithPrefactor[1,(Transpose[Append[tomatrix /@ (\[Alpha]vecs[[k]]),idmatrix],{3,2,1}])],{k,Length[\[Alpha]vecs]}]];

linearconstraints=Table[PositiveMatrixWithPrefactor[1,{{Append[\[Beta],0]}}],{\[Beta],\[Beta]vecs}];

norm=almostzeroarray[Length[\[Alpha]vecs[[1]]]+1,-1->1];
obj=\[DoubleStruckCapitalC]0[Nd,Md]~Join~{0};

WriteBootstrapSDP[fname,SDP[obj,norm,quadraticconstraints~Join~linearconstraints]];];
