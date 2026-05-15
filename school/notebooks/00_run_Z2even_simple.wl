(* ::Package:: *)

(* ::Title:: *)
(*00_run_Z2even_simple: laptop runner for Z2even_Bootstrap_4D_simple.wl*)


(* This file is intentionally a runner, not the bootstrap engine.
   Keep the definitions in Z2even_Bootstrap_4D_simple.wl and use this file to:
   1. load the library;
   2. choose a truncation and options;
   3. run sanity checks;
   4. export a PMP JSON file;
   5. print, or optionally run, the pmp2sdp and sdpb commands.
*)


(* ::Section:: *)
(*1. Locate folder and load the library*)


ClearAll[runnerDir, engineCandidates, engineFile];

runnerDir = Quiet @ Check[NotebookDirectory[], Directory[]];
If[! StringQ[runnerDir], runnerDir = Directory[]];
SetDirectory[runnerDir];

(* Prefer the threshold-option library if it is present; otherwise use the simple library. *)
engineCandidates = FileNameJoin[{runnerDir, #}] & /@ {
   "Z2even_Bootstrap_4D_simple_threshold_option.wl",
   "Z2even_Bootstrap_4D_simple.wl"
   };

engineFile = SelectFirst[engineCandidates, FileExistsQ, Missing["NotFound"]];

If[MatchQ[engineFile, _Missing],
  Print["ERROR: cannot find a bootstrap library. Looked for:"];
  Print /@ engineCandidates;
  Abort[];
];

(* Load definitions only.  Some development versions of the library still contain
   old cluster "Running Functions" and a "Main Section" at the bottom.  We strip
   those here so this runner stays in control. *)
ClearAll[engineSource];
engineSource = Import[engineFile, "Text"];
engineSource = First @ StringSplit[engineSource, "(*Running Functions*)", 2];
engineSource = StringReplace[engineSource, {
   "checkAnsatzDimensions[10,2]" -> "Null",
   "checkConstantPartialWave[2]//N" -> "Null"
   }];

ToExpression[engineSource];

Print["Loaded bootstrap library definitions: ", engineFile];
Print["Working directory: ", runnerDir];
Print["Integral directory: ", integralsdir];
Print["Auxiliary directory: ", auxiliarydir];
If[ValueQ[sdpbPackageFile], Print["SDPB Mathematica package: ", sdpbPackageFile]];



(* ::Section:: *)
(*2. Choose bootstrap options*)


(* Toggle this after the library is loaded.  This only has an effect if the loaded
   library contains the threshold-bound-state option. *)
ClearAll[includeThresholdBoundState];
includeThresholdBoundState = True;

$IncludeThresholdBoundState = includeThresholdBoundState;

Print["Include threshold bound-state basis? ", TrueQ[$IncludeThresholdBoundState]];

(* For a laptop smoke test, keep these small. Increase only after the pipeline works. *)
ClearAll[Nmax, Lmax, theta];
Nmax = 2;
Lmax = 2;
theta = 0;

Print["Using truncation: ", <|"Nmax" -> Nmax, "Lmax" -> Lmax, "theta" -> theta|>];
Print["Number of s-grid points: ", Length[sgrid]];



(* ::Section:: *)
(*3. Positivity t-grid*)


(* The production t-grid in the library can be large.  For the first laptop run,
   use a small grid.  Set useTinyPositivityGrid=False for the production grid. *)
ClearAll[useTinyPositivityGrid, tinyTGrid];
useTinyPositivityGrid = True;
tinyTGrid = SetPrecision[{1/2, 1, 3/2, 2, 5/2, 3, 7/2}, precgrid];

If[TrueQ[useTinyPositivityGrid],
  \[ScriptCapitalT]list = tinyTGrid;
  Print["Using tiny positivity t-grid for smoke test."],
  Print["Using positivity t-grid from the library."]
];

Print["Number of positivity t-points: ", Length[\[ScriptCapitalT]list]];



(* ::Section:: *)
(*4. Check that required integrals exist*)


(*ClearAll[neededIntegralSpecs, integralAvailableQ, missingIntegralSpecs];

neededIntegralSpecs[Nd_Integer, ellMax_Integer] := DeleteDuplicates @ Join[
   Flatten[
    Table[{ell, \[Sigma]center0, n}, {ell, 0, ellMax, 2}, {n, 0, Nd}],
    1
    ],
   Flatten[
    Table[
     If[n + m <= Nd, {ell, \[Sigma]center0, n, m}, Nothing],
     {ell, 0, ellMax, 2}, {n, 1, Nd}, {m, 1, n}
     ],
    2
    ]
   ];

integralAvailableQ[spec_List] := Module[{suffix, hits},
  suffix = integralSuffix @@ spec;
  hits = Select[integralFilesInDirectory[], StringEndsQ[FileBaseName[#], suffix] &];
  hits =!= {}
  ];

missingIntegralSpecs[Nd_Integer, ellMax_Integer] :=
 Select[neededIntegralSpecs[Nd, ellMax], ! integralAvailableQ[#] &];

(* The constant partial-wave check below uses spins 0,2,4, so check up to Max[Lmax,4]. *)
ClearAll[integralCheckLmax, missing];
integralCheckLmax = Max[Lmax, 4];
missing = missingIntegralSpecs[Nmax, integralCheckLmax];

If[missing === {},
  Print["All required integral files are present up to spin ", integralCheckLmax, "."],
  Print["Missing ", Length[missing], " required integral files. First few missing specs:"];
  Print[Take[missing, UpTo[20]]];
  Print["Integral directory: ", integralsdir];
  Abort[];
];

(* Very basic raw integral sanity check. *)
ClearAll[raw0];
raw0 = preintReg[0][\[Sigma]center0, 0];
Print["Raw integral check {Length[raw0], Length[sgrid], MaxAbs[raw0-2]} = ",
  {Length[raw0], Length[sgrid], N[Max[Abs[raw0 - 2]]]}
];
*)


ClearAll[
  neededIntegralSpecs, integralAvailableQ, missingIntegralSpecs,
  downloadAndInstallIntegralCache
];

$IntegralCacheURL =
  "https://github.com/andreaguerr/Smatrix_4d_primal_gapped_singlet/releases/download/integrals-refined-4d-v2026-05-15/integral_storage_bin__refined-4d-gapped__v2026-05-15.tar.gz";

$IntegralCacheArchiveName =
  "integral_storage_bin__refined-4d-gapped__v2026-05-15.tar.gz";

neededIntegralSpecs[Nd_Integer, ellMax_Integer] := DeleteDuplicates @ Join[
   Flatten[
    Table[{ell, \[Sigma]center0, n}, {ell, 0, ellMax, 2}, {n, 0, Nd}],
    1
    ],
   Flatten[
    Table[
     If[n + m <= Nd, {ell, \[Sigma]center0, n, m}, Nothing],
     {ell, 0, ellMax, 2}, {n, 1, Nd}, {m, 1, n}
     ],
    2
    ]
   ];

integralAvailableQ[spec_List] := Module[{suffix, hits},
  suffix = integralSuffix @@ spec;
  hits = Select[integralFilesInDirectory[], StringEndsQ[FileBaseName[#], suffix] &];
  hits =!= {}
  ];

missingIntegralSpecs[Nd_Integer, ellMax_Integer] :=
 Select[neededIntegralSpecs[Nd, ellMax], ! integralAvailableQ[#] &];

downloadAndInstallIntegralCache[] := Module[
  {parentDir, archiveFile, tmpDir, result, sourceDirs, sourceDir, files},

  If[! ValueQ[integralsdir],
   Print["integralsdir is not defined."];
   Abort[];
   ];

  parentDir = DirectoryName[integralsdir];

  If[! DirectoryQ[parentDir],
   CreateDirectory[parentDir, CreateIntermediateDirectories -> True]
   ];

  If[! DirectoryQ[integralsdir],
   CreateDirectory[integralsdir, CreateIntermediateDirectories -> True]
   ];

  archiveFile = FileNameJoin[{parentDir, $IntegralCacheArchiveName}];

  If[! FileExistsQ[archiveFile],
   Print["Downloading integral archive:"];
   Print[$IntegralCacheURL];
   Print["to"];
   Print[archiveFile];

   result = Quiet @ Check[
      URLDownload[$IntegralCacheURL, archiveFile],
      $Failed
      ];

   If[result === $Failed || ! FileExistsQ[archiveFile],
    Print["URLDownload failed. Trying curl..."];
    result = RunProcess[
      {"curl", "-L", $IntegralCacheURL, "-o", archiveFile},
      "ExitCode"
      ];

    If[result =!= 0 || ! FileExistsQ[archiveFile],
     Print["Failed to download integral archive."];
     Abort[];
     ];
    ];
   ,
   Print["Integral archive already downloaded: ", archiveFile];
   ];

  tmpDir = CreateDirectory[];

  Print["Extracting archive..."];

  result = Quiet @ Check[
     ExtractArchive[archiveFile, tmpDir],
     $Failed
     ];

  If[result === $Failed,
   Print["ExtractArchive failed. Trying tar..."];
   result = RunProcess[
     {"tar", "-xzf", archiveFile, "-C", tmpDir},
     "ExitCode"
     ];
   
   If[result =!= 0,
    Print["Failed to extract integral archive."];
     Abort[];
     ];
   ];

  sourceDirs =
   DeleteDuplicates[
    DirectoryName /@ FileNames["rhorhoint*.bin", tmpDir, Infinity]
    ];

  If[sourceDirs === {},
   Print["Could not find rhorhoint*.bin files in extracted archive."];
   Abort[];
   ];

  sourceDir = First @ MaximalBy[
     sourceDirs,
     Length[FileNames["rhorhoint*.bin", #]] &
     ];

  files = FileNames["rhorhoint*.bin", sourceDir];

  Print["Installing ", Length[files], " integral files into:"];
  Print[integralsdir];

  Do[
   CopyFile[
    file,
    FileNameJoin[{integralsdir, FileNameTake[file]}],
    OverwriteTarget -> True
    ],
   {file, files}
   ];

  Print["Integral installation complete."];
  ];


(* The constant partial-wave check below uses spins 0,2,4, so check up to Max[Lmax,4]. *)

ClearAll[integralCheckLmax, missing];

integralCheckLmax = Max[Lmax, 4];
missing = missingIntegralSpecs[Nmax, integralCheckLmax];

If[missing === {},
  Print["All required integral files are present up to spin ", integralCheckLmax, "."],
  
  Print["Missing ", Length[missing], " required integral files. First few missing specs:"];
  Print[Take[missing, UpTo[20]]];
  Print["Integral directory: ", integralsdir];

  downloadAndInstallIntegralCache[];

  missing = missingIntegralSpecs[Nmax, integralCheckLmax];

  If[missing === {},
   Print["All required integral files are now present up to spin ", integralCheckLmax, "."],

   Print["Still missing ", Length[missing], " required integral files after download. First few:"];
   Print[Take[missing, UpTo[20]]];
   Print["Integral directory: ", integralsdir];
   Abort[];
   ];
  ];

(* Very basic raw integral sanity check. *)

ClearAll[raw0];

raw0 = preintReg[0][\[Sigma]center0, 0];

Print[
 "Raw integral check {Length[raw0], Length[sgrid], MaxAbs[raw0-2]} = ",
 {Length[raw0], Length[sgrid], N[Max[Abs[raw0 - 2]]]}
 ];


(* ::Section:: *)
(*5. Run consistency checks*)


Print["Ansatz / partial-wave dimension check:"];
Print[checkAnsatzDimensions[Nmax, Lmax]];

Print["Constant-amplitude partial-wave check:"];
Print[N[checkConstantPartialWave[Nmax]]];

(* Also check the three central objects agree in length. *)
Print["Basis lengths {alpha variables, amplitude basis, Tpw columns}:"];
Print[{Length[\[Alpha]abc[Nmax]], Length[Amplitude[s, t, u][Nmax]], Dimensions[Tpw[0][Nmax]][[2]]}];



(* ::Section:: *)
(*6. Export PMP JSON*)


ClearAll[runsDir, runTag, thresholdTagLocal, baseName, jsonFile];
runsDir = FileNameJoin[{runnerDir, "runs"}];
If[! DirectoryQ[runsDir], CreateDirectory[runsDir, CreateIntermediateDirectories -> True]];

thresholdTagLocal = If[TrueQ[$IncludeThresholdBoundState], "withThreshold", "noThreshold"];
runTag = "Z2evenSimple_" <> thresholdTagLocal <>
  "_N" <> ToString[Nmax] <>
  "_L" <> ToString[Lmax] <>
  "_theta" <> StringReplace[ToString[InputForm[N[theta, 8]]], {"." -> "p", "-" -> "m", "`" -> ""}];

baseName = FileNameJoin[{runsDir, runTag}];

jsonFile = buildSmatrixProgram[baseName, Nmax, Lmax, theta];

Print["Wrote PMP JSON: ", jsonFile];
Print["File exists? ", FileExistsQ[jsonFile]];



(* ::Section:: *)
(*7. SDPB options handle and command generation*)


(* This is the only place where SDPB runtime options should be edited.
   DecimalPrecision is in decimal digits, as used by the Mathematica export.
   SDPB itself wants binary bits; we convert below. *)

ClearAll[$SDPBOptions];

$SDPBOptions = <|
  "DecimalPrecision" -> sdpbPrec,
  "UseMPI" -> True,
  "Cores" -> 4,
  "Pmp2SdpExecutable" -> "pmp2sdp",
  "SDPBExecutable" -> "sdpb",
  "MpirunExecutable" -> "mpirun",
  "CheckpointInterval" -> 50000,
  "MaxRuntime" -> 1072800,
  "MaxIterations" -> 100000,
  "InitialMatrixScalePrimal" -> "1e30",
  "InitialMatrixScaleDual" -> "1e30",
  "DualityGapThreshold" -> "1e-12",
  "MaxComplementarity" -> "1e100",
  "ProcsPerNode" -> Automatic,
  "Environment" -> {
    "OMP_NUM_THREADS=1",
    "OPENBLAS_NUM_THREADS=1"
  },
  "ExtraPmp2SdpOptions" -> {},
  "ExtraSDPBOptions" -> {}
|>;

(* Local machine override. Edit this block only. *)
$SDPBOptions["DecimalPrecision"] = sdpbPrec;
$SDPBOptions["Cores"] = Min[8, $ProcessorCount];
$SDPBOptions["UseMPI"] = True;

(* If binaries are not in your PATH, uncomment and edit:
$SDPBOptions["Pmp2SdpExecutable"] = "/path/to/pmp2sdp";
$SDPBOptions["SDPBExecutable"] = "/path/to/sdpb";
*)

ClearAll[shellQuote, sdpbPrecisionBits];

shellQuote[str_String] := "'" <> StringReplace[str, "'" -> "'\\''"] <> "'";

sdpbPrecisionBits[opts_: $SDPBOptions] :=
  Ceiling[opts["DecimalPrecision"] Log[2, 10]];

ClearAll[pmp2sdpCommand];

pmp2sdpCommand[jsonFile_String, sdpDir_String, opts_: $SDPBOptions] :=
 Module[{prec, args},
  prec = sdpbPrecisionBits[opts];
  args = Join[
    {
     opts["Pmp2SdpExecutable"],
     "--precision=" <> ToString[prec],
     "--input=" <> shellQuote[jsonFile],
     "--output=" <> shellQuote[sdpDir]
    },
    opts["ExtraPmp2SdpOptions"]
  ];
  StringRiffle[args, " "]
 ];

ClearAll[sdpbCommand];

sdpbCommand[sdpDir_String, outDir_String, ckDir_String, opts_: $SDPBOptions] :=
 Module[{prec, envPrefix, mpiPrefix, procsPerNode, args},
  prec = sdpbPrecisionBits[opts];

  envPrefix =
   If[Length[opts["Environment"]] > 0,
    StringRiffle[opts["Environment"], " "] <> " ",
    ""
   ];

  mpiPrefix =
   If[TrueQ[opts["UseMPI"]],
    {opts["MpirunExecutable"], "-n", ToString[opts["Cores"]]},
    {}
   ];

  procsPerNode =
   If[opts["ProcsPerNode"] === Automatic,
    {},
    {"--procsPerNode", ToString[opts["ProcsPerNode"]]}
   ];

  args = Join[
    mpiPrefix,
    {
     opts["SDPBExecutable"],
     "--precision", ToString[prec],
     "-s", shellQuote[sdpDir],
     "-o", shellQuote[outDir],
     "-c", shellQuote[ckDir],
     "--checkpointInterval", ToString[opts["CheckpointInterval"]],
     "--maxRuntime", ToString[opts["MaxRuntime"]],
     "--maxIterations", ToString[opts["MaxIterations"]],
     "--initialMatrixScalePrimal", ToString[opts["InitialMatrixScalePrimal"]],
     "--initialMatrixScaleDual", ToString[opts["InitialMatrixScaleDual"]],
     "--dualityGapThreshold", ToString[opts["DualityGapThreshold"]],
     "--maxComplementarity", ToString[opts["MaxComplementarity"]]
    },
    procsPerNode,
    opts["ExtraSDPBOptions"]
  ];

  envPrefix <> StringRiffle[args, " "]
 ];

ClearAll[sdpDir, outDir, ckDir, logFile, pmp2sdpCmd, sdpbCmd];

sdpDir = StringReplace[jsonFile, ".json" ~~ EndOfString -> "_sdp"];
outDir = StringReplace[jsonFile, ".json" ~~ EndOfString -> "_out"];
ckDir = StringReplace[jsonFile, ".json" ~~ EndOfString -> "_ck"];
logFile = StringReplace[jsonFile, ".json" ~~ EndOfString -> ".log"];

pmp2sdpCmd = pmp2sdpCommand[jsonFile, sdpDir];
sdpbCmd = sdpbCommand[sdpDir, outDir, ckDir];

Print["SDPB options:"];
Print[$SDPBOptions];
Print["SDPB binary precision bits: ", sdpbPrecisionBits[]];

Print["Run this in a terminal to convert JSON to SDPB input:"];
Print[pmp2sdpCmd];
Print["Then run SDPB:"];
Print[sdpbCmd <> " > " <> shellQuote[logFile] <> " 2>&1"];


(* ::Section:: *)
(*8. Optional: run SDPB from Mathematica*)


ClearAll[runPmp2SdpFromMathematica, runSdpbFromMathematica];
runPmp2SdpFromMathematica = False;
runSdpbFromMathematica = False;

If[TrueQ[runPmp2SdpFromMathematica],
  Print["Running pmp2sdp..."];
  Print[RunProcess[{"bash", "-lc", pmp2sdpCmd}, "StandardOutput"]];
];

If[TrueQ[runSdpbFromMathematica],
  Print["Running sdpb..."];
  Print[RunProcess[{"bash", "-lc", sdpbCmd}, "StandardOutput"]];
];



(* ::Section:: *)
(*9. Next steps*)


Print["Smoke-test suggestion: first run Nmax=2, Lmax=0. Then try Lmax=4. Then try Nmax=4."];
Print["To include the threshold bound-state basis, set includeThresholdBoundState=True near the top and rerun from the beginning."];

