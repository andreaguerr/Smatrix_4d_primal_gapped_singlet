(* ::Package:: *)

(* ::Section:: *)
(*Reading Parameters*)


(* ::Subsection::Closed:: *)
(*Basic constants and paths*)


ClearAll[dir];
dir = SetDirectory[NotebookDirectory[]];

ClearAll[precgrid, prec, sdpbPrec];
precgrid = 350;
prec = 250;
sdpbPrec = 250;
$MaxExtraPrecision = 100000;

ClearAll[envPath, firstExistingDirectory];
envPath[name_String] := Module[{v = Environment[name]}, If[StringQ[v] && v =!= "", v, Missing["NotSet"]]];
firstExistingDirectory[candidates_List, fallback_String] := Module[{hits},
  hits = Select[candidates, StringQ[#] && DirectoryQ[#] &];
  If[hits === {}, fallback, First[hits]]
];

ClearAll[integralsdir, auxiliarydir, AllowIntegralGeneration, IntegralFilePrefixes];
integralsdir = firstExistingDirectory[
   DeleteMissing[{envPath["SMATRIX_INTEGRALS_DIR"], 
   FileNameJoin[{dir, "integral_storage_bin"}]}],
   FileNameJoin[{dir, "integral_storage_bin"}]
];
auxiliarydir = firstExistingDirectory[
   DeleteMissing[{envPath["SMATRIX_AUXILIARY_DIR"], 
   FileNameJoin[{dir, "auxiliary_storage_bin"}]}],
   FileNameJoin[{dir, "auxiliary_storage_bin"}]
];

If[! DirectoryQ[integralsdir], CreateDirectory[integralsdir, CreateIntermediateDirectories -> True]];
If[! DirectoryQ[auxiliarydir], CreateDirectory[auxiliarydir, CreateIntermediateDirectories -> True]];
If[! DirectoryQ[FileNameJoin[{auxiliarydir, "positivity_constraints"}]], CreateDirectory[FileNameJoin[{auxiliarydir, "positivity_constraints"}], CreateIntermediateDirectories -> True]];

AllowIntegralGeneration = False;
IntegralFilePrefixes = {"rhorhoint_refined_4d_gapped", "rhorhoint_4d_gapped"};

(* Load official SDPB Mathematica package. Keep SDPB.m next to this file. *)
ClearAll[sdpbPackageFile];
sdpbPackageFile = FileNameJoin[{dir, "SDPB.m"}];
If[FileExistsQ[sdpbPackageFile],
  Get[sdpbPackageFile],
  Print["WARNING: SDPB.m not found at ", sdpbPackageFile, ". Download it from the official SDPB repository before exporting."]
];



(* ::Section:: *)
(*Setup (this code generates the input file)*)


(* ::Subsection::Closed:: *)
(*Coordinates and ansatz variables*)


ClearAll[t,u,\[Rho]\[Rho],\[ScriptS]];

t[s_,x_]:=1/2 (-4+s) (-1+x); 
u[s_,x_]:=-(1/2) (-4+s) (1+x);

\[Rho]\[Rho][s_,\[Sigma]_]:=(Sqrt[\[Sigma]-4]-Sqrt[4-s])/(Sqrt[\[Sigma]-4]+Sqrt[4-s]);
\[ScriptS][\[Rho]_,\[Sigma]_]:=-((-8+\[Rho]^2 (-8+\[Sigma])+\[Sigma]-2 \[Rho] \[Sigma])/(1+\[Rho])^2);


ClearAll[\[Sigma]center0]
\[Sigma]center0=20/3;


(* Cache/grid metadata: must match the downloaded integral cache *)

ClearAll[\[Sigma]center0, \[Sigma]centerlist, lengthcenters];

\[Sigma]center0 = 20/3;

(* Keep this because the cache filenames and sgrid were built with it *)
\[Sigma]centerlist = {30, 50, 86};
lengthcenters = Length[\[Sigma]centerlist];

ClearAll[npts0, npts1, npts, \[Phi]gridCheb0, \[Phi]gridCheb1, \[Rho]gridCheb0, \[Rho]gridCheb1, sgrid];

npts0 = 300;
npts1 = 150;
npts = npts0 + npts1 lengthcenters;

\[Phi]gridCheb0 =SetPrecision[\[Pi]/2 (1 + Cos[\[Pi] Range[npts0 + 1]/(npts0 + 1)]) // Drop[#, -1] &,precgrid] // Reverse;
\[Phi]gridCheb1 =SetPrecision[\[Pi]/2 (1 + Cos[\[Pi] Range[npts1 + 1]/(npts1 + 1)]) // Drop[#, -1] &,precgrid] // Reverse;

\[Rho]gridCheb0 = SetPrecision[Exp[I \[Phi]gridCheb0], precgrid];
\[Rho]gridCheb1 = SetPrecision[Exp[I \[Phi]gridCheb1], precgrid];

sgrid =Join[Chop[\[ScriptS][\[Rho]gridCheb0, \[Sigma]center0], 10^-prec],Flatten[Table[Chop[\[ScriptS][\[Rho]gridCheb1, \[Sigma]], 10^-prec],{\[Sigma], \[Sigma]centerlist}],1]] // Sort;
(*Length[sgrid]*)


(* ::Subsection:: *)
(*Amplitude*)


ClearAll[\[Rho]disc,\[Rho]ddisc];
\[Rho]disc[s_,t_,u_][\[Sigma]_,n_]:=\[Rho]\[Rho][s,\[Sigma]]^n+\[Rho]\[Rho][t,\[Sigma]]^n+\[Rho]\[Rho][u,\[Sigma]]^n;
\[Rho]ddisc[s_,t_,u_][\[Sigma]_,n_,m_]:=\[Rho]\[Rho][s,\[Sigma]]^n (\[Rho]\[Rho][t,\[Sigma]]^m+\[Rho]\[Rho][u,\[Sigma]]^m)+\[Rho]\[Rho][s,\[Sigma]]^m (\[Rho]\[Rho][t,\[Sigma]]^n+\[Rho]\[Rho][u,\[Sigma]]^n)+(\[Rho]\[Rho][u,\[Sigma]]^n \[Rho]\[Rho][t,\[Sigma]]^m+\[Rho]\[Rho][u,\[Sigma]]^m \[Rho]\[Rho][t,\[Sigma]]^n);


ClearAll[$IncludeThresholdBoundState];
$IncludeThresholdBoundState = True;

ClearAll[ThresholdBasis];
ThresholdBasis[s_, t_, u_] :=1/(\[Rho]\[Rho][s, \[Sigma]center0] - 1)+1/(\[Rho]\[Rho][t, \[Sigma]center0] - 1)+1/(\[Rho]\[Rho][u, \[Sigma]center0] - 1);

ClearAll[\[Alpha]abc];
\[Alpha]abc[Nmax_] := \[Alpha]abc[Nmax, TrueQ[$IncludeThresholdBoundState]];
\[Alpha]abc[Nmax_, includeThreshold_] :=
 Join[
  If[TrueQ[includeThreshold], {th}, {}],
  Table[\[Alpha][n], {n, 0, Nmax}],
  Flatten[
   Table[
    If[n + m <= Nmax, \[Alpha][n, m], Nothing],
    {n, 1, Nmax}, {m, 1, n}
    ],
   1
   ]
  ];


ClearAll[Amplitude];
Amplitude[s_, t_, u_][Nmax_] := Amplitude[s, t, u][Nmax, TrueQ[$IncludeThresholdBoundState]];
Amplitude[s_, t_, u_][Nmax_, includeThreshold_] :=
 Join[
  If[TrueQ[includeThreshold], {ThresholdBasis[s, t, u]}, {}],
  Table[\[Rho]disc[s, t, u][\[Sigma]center0, n], {n, 0, Nmax}],
  Flatten[
   Table[
    If[n + m <= Nmax, \[Rho]ddisc[s, t, u][\[Sigma]center0, n, m], Nothing],
    {n, 1, Nmax}, {m, 1, n}
    ],
   1
   ]
  ];


(* ::Subsection:: *)
(*Partial Wave Integrals*)


ClearAll[integralFilesInDirectory];

integralFilesInDirectory[] :=
 Join[
  FileNames["rhorhoint*.bin", integralsdir],
  FileNames["rhorhoint*.txt", integralsdir]
 ];

ClearAll[integralSuffix];

integralSuffix[l_, \[Sigma]_, inds__Integer] :=
 "_" <> StringRiffle[ToString /@ {l, N[\[Sigma], 4], inds}, "_"];

ClearAll[findIntegralFile];

findIntegralFile[l_, \[Sigma]_, inds__Integer] :=
 Module[{suffix, hits},
  
  suffix = integralSuffix[l, \[Sigma], inds];
  
  hits =
   Select[
    integralFilesInDirectory[],
    StringEndsQ[FileBaseName[#], suffix] &
   ];
  
  If[hits === {},
   Print["Missing precomputed integral: ", {l, \[Sigma], inds}];
   Print["Expected suffix: ", suffix];
   Abort[],
   First[hits]
  ]
 ];


ClearAll[readIntegralFile];

readIntegralFile[file_String] :=
 Which[
  StringEndsQ[file, ".bin"],
  Uncompress[Import[file, "String"]],
  
  StringEndsQ[file, ".txt"],
  ToExpression[Import[file]],
  
  True,
  Print["Unknown integral file extension: ", file];
  Abort[]
 ];


ClearAll[preintReg];

preintReg[l_][\[Sigma]_, n_] :=preintReg[l][\[Sigma], n] =readIntegralFile[findIntegralFile[l, \[Sigma], n]];
preintReg[l_][\[Sigma]_, n_, m_] /; n >= m :=preintReg[l][\[Sigma], n, m] =readIntegralFile[findIntegralFile[l, \[Sigma], n, m]];
preintReg[l_][\[Sigma]_, n_, m_] /; n < m :=preintReg[l][\[Sigma], n, m] =(-1)^l preintReg[l][\[Sigma], m, n];


(* ::Subsection:: *)
(*Partial Waves*)


ClearAll[intDReg, intDDReg];
intDReg[l_][\[Sigma]_, n_] :=1/(32 \[Pi]) (Conjugate[\[Rho]\[Rho][sgrid, \[Sigma]]]^n 2 KroneckerDelta[l, 0]+ preintReg[l][\[Sigma], n]+ (-1)^l preintReg[l][\[Sigma], n]);

intDDReg[l_][\[Sigma]_, n_, m_] :=1/(32 \[Pi]) (Conjugate[\[Rho]\[Rho][sgrid, \[Sigma]]]^n (preintReg[l][\[Sigma], m] + (-1)^l preintReg[l][\[Sigma], m])+ Conjugate[\[Rho]\[Rho][sgrid, \[Sigma]]]^m (preintReg[l][\[Sigma], n] + (-1)^l preintReg[l][\[Sigma], n])+ preintReg[l][\[Sigma], n, m]+ (-1)^l preintReg[l][\[Sigma], n, m]);


ClearAll[TpwTh, Tpw];

TpwTh[\[ScriptL]_] :=
 1/(32 \[Pi]) (
   -2 KroneckerDelta[0, \[ScriptL]]
   + (2 KroneckerDelta[0, \[ScriptL]])/
      (-1 + Conjugate[\[Rho]\[Rho][sgrid, \[Sigma]center0]])
   - (8 Sqrt[2/3]
       (-2 + Sqrt[sgrid])^\[ScriptL]
       (2 + Sqrt[sgrid])^(-1 - \[ScriptL]))/
      (1 + 2 \[ScriptL])
   );

Tpw[\[ScriptL]_][Nmax_] := Tpw[\[ScriptL]][Nmax, TrueQ[$IncludeThresholdBoundState]];
Tpw[\[ScriptL]_][Nmax_, includeThreshold_] :=
 Tpw[\[ScriptL]][Nmax, includeThreshold] =
  Join[
    If[TrueQ[includeThreshold], {TpwTh[\[ScriptL]]}, {}],
    Table[intDReg[\[ScriptL]][\[Sigma]center0, n], {n, 0, Nmax}],
    Flatten[
     Table[
      If[n + m <= Nmax,
       intDDReg[\[ScriptL]][\[Sigma]center0, n, m],
       Nothing
       ],
      {n, 1, Nmax}, {m, 1, n}
      ],
     1
     ]
    ] /. {} -> Nothing // Transpose // Chop[#, 10^-prec] &;


(*Transpose[{sgrid,TpwTh[2]}][[200;;200]]//N*)


(*1/(32 \[Pi])NIntegrate[Amplitude[s,t[s,x],u[s,x]][4]/.s->sgrid[[200]]+\[ImaginaryI] 10^-100,{x,-1,1}]
Tpw[0][4][[200]]//N
%%-%*)


(* ::Subsection:: *)
(*Improved Positivity Constraints [we are imposing "improved" positivity constraints for 0<t<4]*)


ClearAll[ttTag];

ttTag[tt_] :=
 StringReplace[
  ToString[InputForm[N[tt, 16]]],
  {"." -> "p", "-" -> "m", "+" -> "p", " " -> "", "`" -> ""}
  ];


ClearAll[thresholdTag, positivityListFile];

thresholdTag[includeThreshold_: TrueQ[$IncludeThresholdBoundState]] :=
 If[TrueQ[includeThreshold], "withThreshold", "noThreshold"];

positivityListFile[Nd_, tt_, includeThreshold_: TrueQ[$IncludeThresholdBoundState]] := FileNameJoin[{
   auxiliarydir,
   "positivity_constraints",
   "positivity_list_" <>
    thresholdTag[includeThreshold] <> "_" <>
    ToString[Nd] <> "_" <>
    ttTag[tt] <> ".bin"
   }];


ClearAll[PositivityList];

PositivityList[Nd_][tt_] := PositivityList[Nd, TrueQ[$IncludeThresholdBoundState]][tt];

PositivityList[Nd_, includeThreshold_][tt_] :=
 PositivityList[Nd, includeThreshold][tt] =
  Module[{filepath, data, dir},
   
   filepath = positivityListFile[Nd, tt, includeThreshold];
   dir = DirectoryName[filepath];
   
   If[! DirectoryQ[dir],
    CreateDirectory[dir, CreateIntermediateDirectories -> True]
    ];
   
   If[FileExistsQ[filepath],
    data = Uncompress[Import[filepath, "String"]];    
    If[Length[data] =!= Length[sgrid],
     Print["Cached positivity list has wrong length: ", filepath];
     Abort[];
     ];    
    ,
    data =
     Table[
       Im[Conjugate[Amplitude[s, tt, 4 - s - tt][Nd, includeThreshold]]],
       {s, sgrid}
       ] // Chop[#, 10^-(prec/2)] &;
    
    data = SetPrecision[data, 200];
    Export[filepath, Compress[data], "String"];
    ];
   
   data
   ];


ClearAll[chebTGrid];
chebTGrid[n_Integer] := Reverse[SetPrecision[2 (1 + Cos[\[Pi] Range[n + 1]/(n + 1)]),precgrid]];


ClearAll[\[ScriptCapitalT]list];
\[ScriptCapitalT]list =Sort @ DeleteDuplicates @ Join[chebTGrid[10],Select[chebTGrid[2000][[1 ;; -1 ;; 2]],# > 3 + 999/1000 &]];


ClearAll[partialWaveImprovement];
partialWaveImprovement[Nd_, Lmax_][tt_] := partialWaveImprovement[Nd, Lmax, TrueQ[$IncludeThresholdBoundState]][tt];
partialWaveImprovement[Nd_, Lmax_, includeThreshold_][tt_] :=
 Sum[
  16 \[Pi] (2 \[ScriptL] + 1)LegendreP[\[ScriptL], 1 + 2 tt/(sgrid - 4)] Im[Tpw[\[ScriptL]][Nd, includeThreshold]],{\[ScriptL], 0, Lmax, 2}];


ClearAll[ImprovedPositivityList];
ImprovedPositivityList[Nd_, Lmax_][tt_] := ImprovedPositivityList[Nd, Lmax, TrueQ[$IncludeThresholdBoundState]][tt];
ImprovedPositivityList[Nd_, Lmax_, includeThreshold_][tt_] :=
 ImprovedPositivityList[Nd, Lmax, includeThreshold][tt] =
  Chop[
   PositivityList[Nd, includeThreshold][tt] - partialWaveImprovement[Nd, Lmax, includeThreshold][tt],
   10^-(prec/2)
   ];


(* ::Subsection::Closed:: *)
(*Objective Definition*)


ClearAll[\[DoubleStruckCapitalC]0, \[DoubleStruckCapitalC]2, Objective];
\[DoubleStruckCapitalC]0[Nmax_] := \[DoubleStruckCapitalC]0[Nmax, TrueQ[$IncludeThresholdBoundState]];
\[DoubleStruckCapitalC]0[Nmax_, includeThreshold_] :=
 \[DoubleStruckCapitalC]0[Nmax, includeThreshold] =
  Module[{s, t},
    1/(32 \[Pi]) Amplitude[s, t, 4 - s - t][Nmax, includeThreshold] /. s -> 4/3 /. t -> 4/3 // Chop[#, 10^-(prec/2)] &
    ] // SetPrecision[#, precgrid] &;

\[DoubleStruckCapitalC]2[Nmax_] := \[DoubleStruckCapitalC]2[Nmax, TrueQ[$IncludeThresholdBoundState]];
\[DoubleStruckCapitalC]2[Nmax_, includeThreshold_] :=
 \[DoubleStruckCapitalC]2[Nmax, includeThreshold] =
  Module[{s, t},
    1/(32 \[Pi]) (1/4 D[Amplitude[s, t, 4 - s - t][Nmax, includeThreshold], {s, 2}]) /. s -> 4/3 /. t -> 4/3 // Chop[#, 10^-(prec/2)] &
    ] // SetPrecision[#, precgrid] &;

Objective[Nmax_][\[Theta]_:0] := Objective[Nmax, TrueQ[$IncludeThresholdBoundState]][\[Theta]];
Objective[Nmax_, includeThreshold_][\[Theta]_:0] :=
 Cos[\[Theta]] \[DoubleStruckCapitalC]0[Nmax, includeThreshold] + Sin[\[Theta]] \[DoubleStruckCapitalC]2[Nmax, includeThreshold];


(* ::Subsection:: *)
(*Quadratic programming*)


ClearAll[almostzeroarray];
almostzeroarray[len_,repl_] :=ReplacePart[ConstantArray[0,len],repl];


(*ClearAll[buildSmatrixProgram];
buildSmatrixProgram[fname_,Nd_,Md_,Pd_,Lmax_][\[Theta]_]:=
Module[{obj,norm,\[Alpha]vecs,\[Eta]vecs,\[Beta]vecs,\[Beta]vecsUP,\[Beta]vecsDN,quadraticconstraints,linearconstraints,tomatrix,idmatrix},

	tomatrix[a_]:={{Re[a],-Im[a]},{-Im[a],-Re[a]}};
	idmatrix:={{1,1},{1,1}};

	\[Alpha]vecs=Table[Print[\[ScriptL]];
	If[\[ScriptL]<4,
	(Sqrt[(sgrid-4)/sgrid]Tpw[\[ScriptL]][Nd,Md,Pd])[[2 ;; ;; 2]],
	(Sqrt[(sgrid-4)/sgrid]Tpw[\[ScriptL]][Nd,Md,Pd])[[4 ;; ;; 4]]
	],{\[ScriptL],0,Lmax,2}]//Flatten[#,1]&;
	Print[\[Alpha]vecs//Total//Total//N];
	Print[\[Alpha]vecs//Dimensions];

	Print["End quadratic constraints"];
	\[Beta]vecs=Table[Print[tt//N];(ImprovedPositivityList[Nd,Md,Pd,Lmax][tt])[[20 ;; ;; 5]],{tt,\[ScriptCapitalT]list[[;; ;; 2]]}]//Flatten[#,1]&;
	Print[\[Beta]vecs//Total//Total//N];
	Print["End linear constraints"];
	Print[\[Beta]vecs//Dimensions];
	
quadraticconstraints =
	Join[Table[
	PositiveMatrixWithPrefactor[1,(Transpose[Append[tomatrix /@ (\[Alpha]vecs[[k]]),idmatrix],{3,2,1}])],{k,Length[\[Alpha]vecs]}]];

	linearconstraints=Table[PositiveMatrixWithPrefactor[1,{{\[Beta]~Join~{0}}}],{\[Beta],\[Beta]vecs}];

	norm=almostzeroarray[Length[\[Alpha]vecs[[1]]]+1,-1->1];
	(*obj= sign SetPrecision[matsolve[Nd,Md][a00,a11,a20,z0] . T20[Nd,Md][z2],precgrid];*)
	obj= Object[Nd,Md,Pd][\[Theta]]~Join~{0};
	WriteBootstrapSDP[fname,SDP[obj,norm,quadraticconstraints~Join~linearconstraints]];
];*)


ClearAll[almostzeroarray];
almostzeroarray[len_, repl_] := ReplacePart[ConstantArray[0, len], repl];

ClearAll[toJsonName];
toJsonName[fname_String] := If[StringEndsQ[fname, ".json"], fname, fname <> ".json"];

ClearAll[buildSmatrixProgram];
buildSmatrixProgram[fname_, Nd_, Lmax_, \[Theta]_:0] := Module[
  {jsonName, obj, norm, \[Alpha]vecs, \[Beta]vecs, quadraticconstraints, linearconstraints, tomatrix, idmatrix},

  (* Keep the original Single_Component elastic-unitarity convention. Do not replace this by the Z2even matrix unless changing the physics deliberately. *)
  tomatrix[a_] := {{-(1/2) Im[a], Re[a]}, {Re[a], 2 Im[a]}};
  idmatrix := {{1, 0}, {0, 0}};

  \[Alpha]vecs = Table[Sqrt[(sgrid - 4)/sgrid] Tpw[\[ScriptL]][Nd], {\[ScriptL], 0, Lmax, 2}] // Flatten[#, 1] &;
  \[Beta]vecs = ParallelTable[Print[tt // N]; ImprovedPositivityList[Nd, Lmax][tt], {tt, \[ScriptCapitalT]list}] // Flatten[#, 1] &;

  quadraticconstraints = Table[
    PositiveMatrixWithPrefactor[1, Transpose[Append[tomatrix /@ (\[Alpha]vecs[[k]]), idmatrix], {3, 2, 1}]],
    {k, Length[\[Alpha]vecs]}
  ];

  linearconstraints = Table[PositiveMatrixWithPrefactor[1, {{Append[\[Beta], 0]}}], {\[Beta], \[Beta]vecs}];

  norm = almostzeroarray[Length[\[Alpha]vecs[[1]]] + 1, -1 -> 1];
  obj = Objective[Nd][\[Theta]] ~Join~ {0};

  jsonName = toJsonName[fname];
  If[Length[Names["WritePmpJson"]] > 0,
    WritePmpJson[jsonName, SDP[obj, norm, quadraticconstraints ~Join~ linearconstraints], sdpbPrec],
    Print["ERROR: WritePmpJson is not defined. Load the official SDPB.m first."]; Abort[]
  ];
  jsonName
];



(* ::Subsection::Closed:: *)
(*Threshold-state option*)


(*
Set $IncludeThresholdBoundState before evaluating cached objects.
  $IncludeThresholdBoundState = False;  (* default: no threshold basis *)
  $IncludeThresholdBoundState = True;   (* prepend the threshold basis th *)
The memoization keys for Tpw, PositivityList, ImprovedPositivityList and the objective include this switch.
*)



(* ::Subsection::Closed:: *)
(*Sanity checks*)


(*ClearAll[checkAnsatzDimensions, checkConstantPartialWave];
checkAnsatzDimensions[Nd_, Lmax_:0] := Module[{nvar, ntpw},
  nvar = Length[\[Alpha]abc[Nd]];
  ntpw = Dimensions[Tpw[Lmax][Nd]][[2]];
  <|"Variables" -> nvar, "TpwColumns" -> ntpw, "MatchQ" -> (nvar == ntpw)|>
];

ClearAll[checkConstantPartialWave];

checkConstantPartialWave[Nd_] := Module[{vars, col},
  vars = \[Alpha]abc[Nd];
  col = First @ First @ Position[vars, \[Alpha][0]];

  <|
   "constantColumn" -> col,
   "spin0ConstantError" ->
    Max[Abs[Tpw[0][Nd][[All, col]] - 3/(16 \[Pi])]],
   "spin2ConstantLeakage" ->
    Max[Abs[Tpw[2][Nd][[All, col]]]],
   "spin4ConstantLeakage" ->
    Max[Abs[Tpw[4][Nd][[All, col]]]]
   |>
  ];*)


(*checkAnsatzDimensions[10,2]*)


(*checkConstantPartialWave[2]//N*)
