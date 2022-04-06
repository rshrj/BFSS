(* ::Package:: *)

BeginPackage["BFSS`"];


(* Check if the package has already been loaded, in case Get was used instead of Needs. *)

If[ValueQ[BFSS`Private`AlreadyLoaded],
    (* Then *)
    (* Unprotect and clear all symbols, so they can be redefined. Useful
         for debugging, or to reload the package after updating. *)
    Unprotect["BFSS`*"];
    Unprotect["BFSS`Private`*"];
    (* Keep the tensor objects and settings created so far during the
         session, so they don't get deleted when the package is reloaded. *)
    BFSSTemp`MatricesData = BFSS`Private`MatricesData;
    ClearAll["BFSS`*"];
    ClearAll["BFSS`Private`*"];
    BFSS`Private`MatricesData = BFSSTemp`MatricesData;
    Remove["BFSSTemp`*"];
    BFSS`Private`AlreadyLoaded = True
    ,
    (* Else *)
    BFSS`Private`AlreadyLoaded = True;
    (* Initialize the symbol MatricesData, which is used to store the data
         for the matrices objects, as well as user settings. This is done only 
        on first load. *)
    BFSS`Private`MatricesData = Association[];
];


Null[{
	FRestylePlot,
	FC,
	FA,
	FVarsToMatrix,
	FRandomConfig,
	FLagrangian,
	FLagrangianWithSource,
	FSolve,
	FStartupCheckForUpdates,
	FCheckForUpdates
}];


Begin["`Private`"];


BFSSVersion = "v1.0.0 (2022-04-06)";
BFSSURL = "https://raw.githubusercontent.com/rshrj/BFSS/master/BFSS.m";

(* If the global options LocalSymbol has been previously set, then it will have the Head Association. Otherwise it will have the Head LocalSymbol, and we create it now for later use. *)
If[
    Head[LocalSymbol["BFSSGlobalOptions"]] =!= Association,
(* Then *)
    LocalSymbol["BFSSGlobalOptions"] = Association[];
];
(* Load the global options from the persistent storage. *)
BFSSGlobalOptions = LocalSymbol["BFSSGlobalOptions"];
(* Set the "BFSSVersion" key to the current version of the package. *)
BFSSGlobalOptions["BFSSVersion"] = BFSSVersion;
(* If the option "AutoUpdates" has not been set, or is not a boolean value, we set it now to the default value of True. *)
If[
    !KeyExistsQ[BFSSGlobalOptions, "AutoUpdates"] || !BooleanQ[BFSSGlobalOptions["AutoUpdates"]],
(* Then *)
    BFSSGlobalOptions["AutoUpdates"] = True;
];
(* If the option "AllowOverwrite" has not been set, or is not a boolean value, we set it now to the default value of False. *)
If[
    !KeyExistsQ[BFSSGlobalOptions, "AllowOverwrite"] || !BooleanQ[BFSSGlobalOptions["AllowOverwrite"]],
(* Then *)
    BFSSGlobalOptions["AllowOverwrite"] = False;
];
(* If the option "MaxPrecision" has not been set or is not an integer between 10 to 1000 then set it to the default value of 150 *)
If[
    !KeyExistsQ[BFSSGlobalOptions, "MaxPrecision"] || !IntegerQ[BFSSGlobalOptions["MaxPrecision"]] || !Between[BFSSGlobalOptions["MaxPrecision"], {10,1000}],
(* Then *)
    BFSSGlobalOptions["MaxPrecision"] = 150;
];
(* If the option "Dimensions" has not been set or is not a list {Nd, Kd} with 2<=Nd<=25, 2<=Kd<=25 then set it to the default value of {2, 3} *)
If[
    !KeyExistsQ[BFSSGlobalOptions, "Dimensions"] || !ListQ[BFSSGlobalOptions["Dimensions"]] || Dimensions[BFSSGlobalOptions["Dimensions"]]!={2} || !ForAll[BFSSGlobalOptions["Dimensions"], Between[#, {2, 25}]&],
(* Then *)
    BFSSGlobalOptions["Dimensions"] = {2, 3};
];
(* Save the global options to the persistent storage. *)
LocalSymbol["BFSSGlobalOptions"] = BFSSGlobalOptions;


(* Create a nicely-formatted usage message. *)
CreateUsageMessage[f_, msg_String : {}] := Evaluate[f::usage] = ToString[TextCell[Row[Flatten[{List @@ StringReplace[msg, {"`" ~~ (x : Shortest[__]) ~~ "`" :> Style[x, Bold]}]}]]], StandardForm];


Needs["VariationalMethods`"];


(* Utility functions *)

CreateUsageMessage[FRestylePlot, "FRestylePlot[`plot`, `newStyles`] applies new styles `newStyles` to Graphics object `plot`"];
FRestylePlot[plot_Graphics,styles_List,op:OptionsPattern[Graphics]]:=Module[{x=styles},Show[MapAt[#/. {__,ln__Line}:>{Directive@Last[x=RotateLeft@x],ln}&,plot,1],op]];

CreateUsageMessage[FC, "FC[`A`, `B`] calculates the commutator of `A` and `B`"];
FC[A_,B_]:=A . B-B . A;

CreateUsageMessage[FA, "FA[`A`, `B`] calculates the anti-commutator of `A` and `B`"];
FA[A_,B_]:=A . B+B . A;

CreateUsageMessage[FVarsToMatrix, "FVarsToMatrix[`X`, `t`] creates a list of `Kd` `Nd`x`Nd` matrices with the symbol `X` at time `t`"];
FVarsToMatrix[X_,t_]:=Module[{Nd, Kd}, ({Nd, Kd} = BFSSGlobalOptions["Dimensions"];Table[Table[X[k,a,b][t],{a,1,Nd},{b,1,Nd}],{k,1,Kd}])];



(* Main functions *)

CreateUsageMessage[FRandomConfig, "FRandomConfig[] creates a list of `Kd` `Nd`x`Nd` traceless Hermitian matrices with random entries"];
FRandomConfig[]:=Module[{X,XH,Nd,Kd},(
	{Nd, Kd} = BFSSGlobalOptions["Dimensions"];
	
	X=Table[
		Table[
			RandomComplex[WorkingPrecision->BFSSGlobalOptions["MaxPrecision"]],
		{a,1,Nd},{b,1,Nd}],
	{k,1,Kd}];
	
	XH=Table[
		(X[[k]]+X[[k]]\[ConjugateTranspose])/2-1/Nd Tr[(X[[k]]+X[[k]]\[ConjugateTranspose])/2]*IdentityMatrix[Nd],
	{k,1,Kd}];
	
	XH
)];

CreateUsageMessage[FLagrangian, "FLagrangian[`X`, `t`] gives the symbolic BFSS lagrangian (no sources) for matrix symbols `X` at time `t`"];
FLagrangian[X_,t_]:=Module[{Y=FVarsToMatrix[X,t], Nd, Kd},(
	{Nd, Kd} = BFSSGlobalOptions["Dimensions"];
	1/2 Sum[Tr[D[Y,t][[i]] . D[Y,t][[i]]],{i,1,Kd}] + 1/4 Sum[Tr[FC[Y[[i]],Y[[j]]] . FC[Y[[i]],Y[[j]]]],{i,1,Kd},{j,1,Kd}]
)];

CreateUsageMessage[FLagrangianWithSource, "FLagrangianWithSource[`X`, `t`, `Op`, `j`] gives the symbolic BFSS lagrangian for matrix symbols
	`X` at time `t` with source `j[t]` coupled to the operator `Op[X, t]`"];
FLagrangianWithSource[X_,t_,Op_,j_] := FLagrangian[X, t] + Op[X, t] * j[t]

CreateUsageMessage[FSolve, "FSolve[`X0`, `X0dot`, `X`, {`tStart`, `tEnd`}] attempts to solve the BFSS system with an optional source for
	 initial conditions {`X0`, `X0dot`} between `t` = `tStart` to `tEnd`"];
FSolve::ErrorIncorrectDimensions = "`X0` and/or `X0dot` do not have the right dimensions. Check dimensions or use FSetDimensions.";
Options[FSolve] = {
	"Coupling" -> {0&, 0&},
	"WorkingPrecision" -> 20,
	"PrecisionGoal" -> 15,
	"AccuracyGoal" -> 15,
	"MaxSteps" -> 10^6
};
FSolve[X0_, X0dot_, X_, {tStart_, tEnd_}, OptionsPattern[]] := Module[{Nd, Kd, Op = OptionValue["Coupling"][[1]], j = OptionValue["Coupling"][[2]], eqns, initConditions, t, time}, (
	{Nd, Kd} = BFSSGlobalOptions["Dimensions"];
	
	If[Dimensions[X0] != {Kd, Nd, Nd} || Dimensions[X0dot] != {Kd, Nd, Nd},
	(*Then*)
		Message[FSolve::ErrorIncorrectDimensions];
		Abort[];
	];
	
	eqns = EulerEquations[FLagrangianWithSource[X, t, Op, j], FVarsToMatrix[X, t] // Flatten, t] // Simplify;
	
	initConditions = Table[{X[k, a, b][0] == X0[[k]][a, b], X[k, a, b]'[0] == X0dot[[k]][a, b]}, {k, 1, Kd}, {a, 1, Kd}, {b, 1, Kd}] // Flatten;
	
	Monitor[
		NDSolve[{eqns, initConditions} // Flatten, Table[X[k, a, b], {k, 1, Kd}, {a, 1, Kd}, {b, 1, Kd}] // Flatten, {t, tStart, tEnd},
			StepMonitor :> (time = t),
			PrecisionGoal -> OptionValue["PrecisionGoal"],
			AccuracyGoal -> OptionValue["AccuracyGoal"],
			WorkingPrecision -> OptionValue["WorkingPrecision"],
			MaxSteps -> OptionValue["MaxSteps"]] // First,
		ProgressIndicator[time/(tEnd-tStart)]
	]
)];


FStartupCheckForUpdates[] := Module[
    {
        newVersion,
        remoteFile,
        versionLookup
    },
    remoteFile = Quiet[Import[BFSSURL, "Text"]];
    Unprotect[UpdateMessage];
    If[
        remoteFile === $Failed,
    (* Then *)
        UpdateMessage = Row[{"Could not check for updates automatically. Please visit ", Hyperlink["https://github.com/rshrj/BFSS"], " to check manually."}],
    (* Else *)
        versionLookup = StringCases[remoteFile, Shortest["BFSSVersion = \"" ~~ __ ~~ "\";"]];
        If[
            Length[versionLookup] == 1,
        (* Then *)
            newVersion = StringTake[versionLookup[[1]], {16, StringLength[versionLookup[[1]]] - 2}];
            If[
                newVersion === BFSSVersion,
            (* Then *)
                UpdateMessage = "You have the latest version of the package.",
            (* Else *)
                UpdateMessage = Row[{"A new version of the package is available: ", Style[newVersion, Bold], ". For more information, type ", CreateButton["FCheckForUpdates[]", FCheckForUpdates[]], "."}];
            ];
        ];
    ];
    Protect[UpdateMessage];
];


(* Check for updates at startup, if the AutoUpdates setting is turned on. *)
If[
    BFSSGlobalOptions["AutoUpdates"],
(* Then *)
    UpdateMessage = "Checking for updates...";
    UpdateCheck = Row[{Dynamic[UpdateMessage], " To disable automatic checks for updates at startup, type ", CreateButton["FSetAutoUpdates[False]", FSetAutoUpdates[False]], "."}];
    SessionSubmit[StartupCheckForUpdates[]],
(* Else *)
    UpdateCheck = Row[{"To check for updates, type ", CreateButton["FCheckForUpdates[]", FCheckForUpdates[]], ". ", "To enable automatic checks for updates at startup, type ", CreateButton["FSetAutoUpdates[True]", FSetAutoUpdates[True]], "."}];
]


CreateUsageMessage[FCheckForUpdates, "FCheckForUpdates[] checks the GitHub repository for new versions of this package. If a new version is available, the user will be given the option to download or install it."];
FCheckForUpdates[] := Module[
    {
        errorMessage,
        newVersion,
        remoteFile,
        toPrint,
        versionLookup
    },
    errorMessage = Row[{"Error: Failed to check for updates. Please visit ", Hyperlink["https://github.com/rshrj/BFSS"], " to check manually."}];
    BFSSPrint["Checking GitHub repository for updates..."];
    remoteFile = Quiet[Import[BFSSURL, "Text"]];
    If[
        remoteFile === $Failed,
    (* Then *)
        BFSSPrint[errorMessage];
        Abort[];
    ];
    versionLookup = StringCases[remoteFile, Shortest["BFSSVersion = \"" ~~ __ ~~ "\";"]];
    If[
        Length[versionLookup] == 1,
    (* Then *)
        newVersion = StringTake[versionLookup[[1]], {16, StringLength[versionLookup[[1]]] - 2}];
        If[
            newVersion === BFSSVersion,
        (* Then *)
            BFSSPrint["You have the latest version of the package."],
        (* Else *)
            toPrint = {Row[{"A new version of the package is available: ", Style[newVersion, Bold]}]};
            AppendTo[toPrint, Row[{"\[Bullet] ", CreateButton[
                "Visit GitHub repository.",
                SystemOpen["https://github.com/rshrj/BFSS"]
            ]}]];
            AppendTo[toPrint, Row[{"\[Bullet] ", CreateButton[
                "Reload new version directly from GitHub without downloading it.",
                Get[BFSSURL]
            ]}]];
            (* If the notebook is an Untitled notebook, meaning it is not an actual file in the file system, then NotebookDirectory[] will return $Failed and issue the error message NotebookDirectory::nosv. Otherwise, show an option to download the notebook to the current notebook directory. *)
            Off[NotebookDirectory::nosv];
            If[
                NotebookDirectory[] =!= $Failed,
            (* Then *)
                AppendTo[toPrint, Row[{"\[Bullet] ", CreateButton[
                    "Download new version to " <> FileNameJoin[{NotebookDirectory[], "BFSS.m"}] <> " and reload the package.",
                    URLDownload[BFSSURL, FileNameJoin[{NotebookDirectory[], "BFSS.m"}]]; BFSSPrint["Downloaded! Reloading..."]; Get[FileNameJoin[{NotebookDirectory[], "BFSS.m"}]]
                ]}]];
            ];
            On[NotebookDirectory::nosv];
            AppendTo[toPrint, Row[{"\[Bullet] ", CreateButton[
                "Install new version to " <> FileNameJoin[{$UserBaseDirectory, "Applications", "BFSS.m"}] <> " and reload the package.",
                URLDownload[BFSSURL, FileNameJoin[{$UserBaseDirectory, "Applications", "BFSS.m"}]]; BFSSPrint["Installed! Reloading..."]; Get[FileNameJoin[{$UserBaseDirectory, "Applications", "BFSS.m"}]]
            ]}]];
            BFSSPrint[Column[toPrint]];
        ],
    (* Else *)
        BFSSPrint[errorMessage];
    ];
];


BFSSPrint[expression_] := CellPrint[ExpressionCell[expression, "Output", Editable -> False, CellLabel -> "BFSS:", CellLabelStyle -> Directive["CellLabel", Smaller, Blue]]];
BFSSPrint[expressions__] := BFSSPrint[Row[{expressions}]];

(* Create a clickable button that looks and behaves like a hyperlink. *)
CreateButton[label_, action_] := Button[
    MouseAppearance[Mouseover[
        Style[label, "Hyperlink"],
        Style[label, "HyperlinkActive"]
    ], "LinkHand"],
    action,
    Appearance -> "Frameless",
    BaseStyle -> "Hyperlink"
];
Attributes[CreateButton] = HoldRest;


(* Print a welcome message at startup. *)
BFSSPrint[Column[{
    Style[Row[{"BFSS: A package to do Matrix model simulations in Mathematica"}], Bold, Larger],
    Style[Row[{"By Rishi Raj (", Hyperlink["rishiraj.1012exp@gmail.com", "mailto:rishiraj.1012exp@gmail.com"], ")"}], Bold],
    Style[Row[{BFSSVersion}], Bold],
    Style[Row[{"GitHub repository: ", Hyperlink["https://github.com/rshrj/BFSS"]}], Bold],
    Row[{"\[Bullet] To view the full documentation for the package, type ", CreateButton["FDocs[]", FDocs[]], "."}],
    Row[{"\[Bullet] To list all available modules, type ", CreateButton["?BFSS`*", BFSSPrint[Information["BFSS`*"]]], "."}],
    Row[{"\[Bullet] To get help on a particular module, type ", Style["?", "Input"], " followed by the module name."}],
    Row[{"\[Bullet] ", UpdateCheck}]
}]];


End[]
EndPackage[];
