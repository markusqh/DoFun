(* ::Package:: *)

(* Mathematica Init File *)

BeginPackage["DoFun`"]

$doFunMainVersion=3;
$doFunSubVersion=0;
$doFunBuiltVersion=1;

$doFunVersion=ToString@$doFunMainVersion~~"."~~ToString@$doFunSubVersion~~"."~~ToString@$doFunBuiltVersion;

(* check for new version *)
(* updated in 3.0.1 to use info from Github *)

getLatestVersionNumbers::usage = "Gets the version number of the latest version on GitHub.";

checkForDoFunUpdates::usage = "checkForDoFunUpdates[] searches online for DoFun updates.
checkForDoFunUpdates[quiet] with quiet=True suppresses warnings if paclet info is not found.";

updateDoFun::usage = "updateDoFun[] updates DoFun to the latest release from GitHub.";
  
getLatestVersionNumbers::pacletinfonotfound="Paclet info could not be found at `1`. Ensure that you have a working network connection.";

checkForDoFunUpdates::usage="checkForDoFunUpdates[] searches online for DoFun updates.
checkForDoFunUpdates[quiet] with quiet=True suppresses warnings if paclet info is not found.";

General::allowinternetuse="You have forbidden Mathematica to access the internet. This function requires internet access.";

Begin["`Private`"];


(* location of GitHub repository *)
doFunRepositoryAddress="https://raw.githubusercontent.com/markusqh/DoFun/master/";


(* get the version number of the latest version on GitHub *)
getLatestVersionNumbers[quiet_]:=Module[{newVersionString,pacletInfoLocation=doFunRepositoryAddress<>"DoFun/PacletInfo.m"},

	If[quiet,
		newVersionString=Quiet[Check[Version/.List@@Import[pacletInfoLocation],""]];,
		newVersionString=Quiet[Check[Version/.List@@Import[pacletInfoLocation],Message[getLatestVersionNumbers::pacletinfonotfound,doFunRepositoryAddress];""],{Import::nffil,ReplaceAll::reps,FetchURL::httperr}];
	];

	Return[If[newVersionString==="",{0,0,0},
		First@StringReplace[newVersionString,mainVersion__~~"."~~version__~~"."~~builtVersion__:>ToExpression/@{mainVersion,version,builtVersion}]
		]
	];
];


(* checks for updates on GitHub; checks version number of master branch, not the latest release *)
checkForDoFunUpdates[quiet_:False]:=Module[{newVersionNumbers,newVersionString},
	If[Not["AllowInternetUse" /. SystemInformation["Network"]],If[Not[quiet],Message[checkForDoFunUpdates::allowinternetuse]];Return[];];
		newVersionNumbers=getLatestVersionNumbers[quiet];
		newVersionString=StringJoin[Riffle[ToString/@newVersionNumbers,"."]
	];

	If[newVersionNumbers=!={0,0,0},
		If[newVersionNumbers[[1]]+newVersionNumbers[[2]]/10>$doFunMainVersion+$doFunSubVersion/10,Print["
DoFun version "<>newVersionString<>" is available. You are currently using version "<> $doFunVersion<>".
You are strongly advised to update as the new version may contain bugfixes. You can do so by evaluating updateDoFun[].
Please be aware that syntax changes may have occured. We recommend to read the changelog before updating.
"];,
			If[newVersionNumbers[[1]]+newVersionNumbers[[2]]/10+newVersionNumbers[[3]]/100>$doFunMainVersion+$doFunSubVersion/10+$doFunBuiltVersion/100,Print["
DoFun version "<>newVersionString<>" is available. You are currently using version "<> $doFunVersion<>".
You can update the DoFun by evaluating updateDoFun[].
"];,
				If[Not[quiet],Print["You are already using the latest version of the DoFun, version "<> $doFunVersion<>"."];
				];
			];
		];
	];
];


(* updates DoFun to latest release on GitHub *)
updateDoFun[]:=Module[{newVersionNumbers},
	If[Not["AllowInternetUse" /. SystemInformation["Network"]],Message[updateDoFun::allowinternetuse];Return[];];
	
	newVersionNumbers=getLatestVersionNumbers[False];
	
	If[newVersionNumbers=!={0,0,0},
		If[newVersionNumbers[[1]]+newVersionNumbers[[2]]/10+newVersionNumbers[[3]]/100>$doFunMainVersion+$doFunSubVersion/10+$doFunBuiltVersion/100,
			Import[doFunRepositoryAddress<>"DoFun/DoFunInstaller.m"];,
			Print["You are already using the latest version of the DoFun, version "<> $doFunVersion<>"."];
		];
	];
];


(*check for updates at startup*)
If["AllowInternetUse" /. SystemInformation["Network"], checkForDoFunUpdates[True];];


End[];

EndPackage[];


Needs/@{"DoFun`DoDSERGE`","DoFun`DoAE`","DoFun`DoFR`"}


If[DoFun`$doFunStartMessage=!=False,
	Print["\nDoFun loaded.
\nVersion "<> $doFunVersion<>
"\nReinhard Alkofer, Jens Braun, Anton K. Cyrol, Markus Q. Huber, Jan M. Pawlowski, Kai Schwenzer, 2008-2019
\nDetails at https://github.com/markusqh/DoFun/."];
];
