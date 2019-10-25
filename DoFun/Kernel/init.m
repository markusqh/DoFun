(* ::Package:: *)

(* Mathematica Init File *)

BeginPackage["DoFun`"]

$doFunMainVersion=3;
$doFunSubVersion=0;
$doFunBuiltVersion=1;

$doFunVersion=ToString@$doFunMainVersion~~"."~~ToString@$doFunSubVersion~~"."~~ToString@$doFunBuiltVersion;

(* check for new version *)
onlineVersion=Quiet[Check[doFunVersion/.ToExpression@Import["http://physik.uni-graz.at/~mqh/DoFun/DoFun_version"],"0.0.0"]];
onlineMainVersion=ToExpression@StringCases[onlineVersion, main__ ~~ "." ~~ sub__~~"."~~built__ :> main][[1]];
onlineSubVersion=ToExpression@StringCases[onlineVersion, main__ ~~ "." ~~ sub__~~"."~~built__ :> sub][[1]];
onlineBuiltVersion=ToExpression@StringCases[onlineVersion, main__ ~~ "." ~~ sub__~~"."~~built__ :> built][[1]];

(* check if override is set *)
override=Quiet[Check[
  ToExpression@Import[FileNameJoin[{$UserBaseDirectory, "Applications", 
      "DoFun_newVersionMessage_override"}]], 0]];

(* messages if new version available in descending order *)
If[override!=1,
If[onlineMainVersion > $doFunMainVersion,
 	overrideNew=ChoiceDialog[
 Hyperlink[
  "There is a major new release (" <> onlineVersion <> 
   ") of DoFun available. You have version " <> $doFunVersion <>". Please go to http://http://physik.uni-graz.at/~mqh/DoFun to download it.", 
  "http://http://physik.uni-graz.at/~mqh/DoFun"], {"OK"->0, "Deactivate this message in the future" -> 
 1}];,

If[onlineMainVersion == $doFunMainVersion && onlineSubVersion > $doFunSubVersion,
overrideNew=ChoiceDialog[
 Hyperlink[
  "There is a new release (" <> onlineVersion <> 
   ") of DoFun available. You have version " <> $doFunVersion <> 
   ". Please go to http://physik.uni-graz.at/~mqh/DoFun to download it.", 
  "http://physik.uni-graz.at/~mqh/DoFun"], {"OK"->0, "Deactivate this message in the future" -> 
 1}];,
 
If[onlineMainVersion == $doFunMainVersion && onlineSubVersion == $doFunSubVersion && onlineBuiltVersion > $doFunBuiltVersion,
overrideNew=ChoiceDialog[
 Hyperlink[
  "There is a minor new release (" <> onlineVersion <> 
   ") of DoFun available. You have version " <> $doFunVersion <> 
   ". Please go to http://physik.uni-graz.at/~mqh/DoFun to download it.", 
  "http://physik.uni-graz.at/~mqh/DoFun"], {"OK"->0, "Deactivate this message in the future" -> 
 1}];   
   ]
  ] 
 ]
]

If[overrideNew==1, Export[FileNameJoin[{$UserBaseDirectory, "Applications", 
     "DoFun_newVersionMessage_override"}], 1, "Text"]]

(* function for resetting the override option *)
resetDoFunVersionCheckOverride[]:=Export[FileNameJoin[{$UserBaseDirectory, "Applications", 
     "DoFun_newVersionMessage_override"}], 0, "Text"]
	
EndPackage[];


Needs/@{"DoFun`DoDSERGE`","DoFun`DoAE`","DoFun`DoFR`"}


If[DoFun`$doFunStartMessage=!=False,
	Print["\nDoFun loaded.
\nVersion "<> $doFunVersion<>
"\nReinhard Alkofer, Jens Braun, Anton K. Cyrol, Markus Q. Huber, Jan M. Pawlowski, Kai Schwenzer, 2008-2019
\nDetails at https://github.com/markusqh/DoFun/."];
];
