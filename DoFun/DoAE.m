(* ::Package:: *)

(* Mathematica Package *)

(* (c) and written by Markus Q. Huber *)

(* published under GNU General Public License v3.0 *)





(* version history *)
(*  0.1 (autumn of 2008): first running version
    0.2 (never released):
    	-) introduced DoDSE`$startMessage to allow suppression of the message when loading the package
   		-) inclusion of RGEs
    2.0.0 (25.2.2011): First public release of DoFun via http://www.tpi.uni-jena.de/qfphysics/homepage/mhub/DoFun/.
    2.0.1 (22.1.2013, never released):
   		-) bugfix in errorhandler getAE::incorrectFormat: correction of first condition that all fields have to be defined as such
   		-) bugfix in getAEDo: fields with two indices of the same type are handled correctly now instead of using the same index names (new function removeDoubleIndices)
		-) bugfix in insertMomenta: ops with only P or V work now (analogously to S)
    2.0.2 (22.1.2015): no changes
    2.0.3 (19.1.2017): no changes
    2.0.4 (5.12.2017): no changes
    3.0.0 (31.7.2019):
    	-) added composite operator CO tofunctionality
    	-) momentum routing for more general diagrams than DSEs and RGEs
*)


BeginPackage["DoFun`DoAE`", { "DoFun`DoDSERGE`", "DoFun`DoFR`"}]

(* variable that gives the version of DoDSE *)
$DoAEVersion="3.0.0";

If[Not@FreeQ[Contexts[],"DoFun`"],DoFun`DoAE`$doAEStartMessage=False];
If[(DoFun`DoAE`$doAEStartMessage=!=False),
	Print["Package DoAE loaded.
\nVersion "<> $DoAEVersion <>
"\nJens Braun, Anton K. Cyrol, Reinhard Alkofer, Markus Q. Huber, Jan. M. Pawlowski, Kai Schwenzer, 2008-2019\n
\nDetails at https://github.com/markusqh/DoFun/."];
];


(* ::Section:: *)
(* usages *)


$loopMomenta::usage="$loopMomenta determines the names given to loop momenta in DoFun.

Default: q1, q2, ...\n
These names are protected and should not be used otherwise.
";

addIndices::usage="addIndices[] gives the list of indices known by DoFun.
addIndices[{iname, dummies}] adds one new index, where iname gives the name of the index and dummies is the list of variable names that should be used.
addIndices[{iname1, dummies1}, {iname2, dummies2}}] adds several indices.

Examples:
addIndices[]
addIndices[{son, {a, b, c, d, e, f, g, h, i, j}}]
";

createDummyListUnique::usage="createDummyListUnique[n, t] creates at least n unique variable names for the index type t.
The option dummyNames[t] sets the variable names for the index type t.
";

explicit::usage="explicit is an option getAE and passed on to V, P, S, CO and dR in the result.
";

getAE::usage="getAE[exp, ls, [opts]] transforms a DSE, RGE or composite operator equation expr  with external legs ls into symbolic form into an algebraic expression.
The external legs are given in the following form:
{{field1, symInd1, mom1, inds1}, {field2, symInd2, mom2, inds2}, ...}.
Each individual list correponds to one external leg of the diagram, where fieldi indicates the field name, symbIndi the index in the symbolic form and momi and indsi the momentum and the indices for the algebraic form .\n

Hint:
A useful option is explicit -> False. With the option explicit -> False, the propagators and vertices are not replaced by their algebraic form but with the generic expressions, which, however, already contain all indices.

Example:
This example takes the so-called sunset diagram from the DSE of a field A. It has only one index adj. For illustration purposes the propagator and the vertices are taken as simple expressions.

Here we show the generic structure in terms of propagators and vertices:
defineFieldsSpecific[{A[momentum, adj]}];
getAE[op[S[{A, i1}, {A, r1}, {A, r2}, {A, s1}], P[{A, r1}, {A, s2}], P[{A, r2}, {A, t2}], P[{A, s1}, {A, u2}], V[{A, i2}, {A, s2}, {A, t2}, {A, u2}]], {{A, i1, p1, a}, {A, i2, p2, b}}, explicit -> False]

And here we replace the generic propagators and vertices by algebraic expressions and simplify the result with integrateDeltas:
defineFieldsSpecific[{A[momentum, adj]}];
P[A[p1_, i1_], A[p2_, i2_], explicit -> True] := delta[adj, i1, i2]/p^2;
S[A[p1_, i1_], A[p2_, i2_], A[p3_, i3_], A[p4_, i4_], explicit -> True] := g delta[adj, i1, i2] delta[adj, i3, i4];
V[A[p1_, i1_], A[p2_, i2_], A[p3_, i3_], A[p4_, i4_], explicit -> True] :=  S[A[p1, i1], A[p2, i2], A[p3, i3], A[p4, i4], explicit -> True];
getAE[op[S[{A, i1}, {A, r1}, {A, r2}, {A, s1}], P[{A, r1}, {A, s2}], P[{A, r2}, {A, t2}], P[{A, s1}, {A, u2}], V[{A, i2}, {A, s2}, {A, t2}, {A, u2}]], {{A, i1, p1, a}, {A, i2, p2, b}}] // integrateDeltas
";

loadFeynCalc::usage="loadFeynCalc[] fixes a problem when loading FeynCalc. Furthermore the output format is set to StandardForm.
loadFeynCalc[pack] with pack the path of the FeynCalc package. If none is given, HighEnergyPhysics`FeynCalc` is used.
This should be used when problems with FeynCalc occur.
";

removeIndices::usage="removeIndices[iname] removes the index iname from the list of known indices.
removeIndices[{iname1, iname2}] removes the indices with names iname1 and iname2.

Example:
resetIndices[]
removeIndices[lor]
";

resetIndices::usage="resetIndices[] resets the known indices to the standard, i.e., only Lorentz and adjoint indices, lor and adj, respectively.
";

save::usage="save is an option of getAE. If set to True, it saves the results which can speed up repeated calculations. Should be used with care.
Default: False.
";




(* ::Section:: *)
(* options and misc *)


Options[getAE]:={save:>False, explicit->True};

Options[createDummyListUnique]={dummyNames[adj]->{Global`a,Global`b,Global`c,Global`d,Global`e,Global`f,Global`g,Global`h,Global`i,Global`j},
	dummyNames[lor]->{Global`\[Mu],Global`\[Nu],Global`\[Rho],Global`\[Sigma],Global`\[Tau],Global`\[Alpha],Global`\[Beta],Global`\[Gamma],Global`\[Delta],Global`\[Epsilon]}};


(* backwards compatibility and abbreviations *)
symb2Alg=getAE;
getAlgebraicExpression=getAE;



Begin["`Private`"]

(* private functions (alphabetic)
	addMomentum
	checkIndicesAlg
	checkMomentumConservation
	contractIndices
	flowDiagram
	flowProps
	getAEDo
	indicesTestAlg
	insertMomenta
	putIndices
	removeDoubleIndices
	replaceIndices

*)


(* names for loop momenta; six should be enough *)

$loopMomenta = Module[{i},Table[ToExpression["Global`q" <> ToString@i], {i, 1, 100}]];
Protect/@$loopMomenta;


(* ::Section:: *)
(* FeynCalc compatibility *)


(* FeynCalc deletes contexts from $ContextPath with the line
Get[Last[HighEnergyPhysics`FeynCalc`Private`tarcerfilenames]];
add the contexts again after loading FeynCalc *)

loadFeynCalc[package_:"HighEnergyPhysics`FeynCalc`"]:=Module[{$contextPathTemp},
	
$contextPathTemp = $ContextPath;
Get@package;

(* add the contexts deleted by FeynCalc *)
$ContextPath = Union[$ContextPath, $contextPathTemp];

(* reset output type *)
SetOptions[$FrontEnd, 
 "CommonDefaultFormatTypes" -> {"Output" -> StandardForm}];
];



(* ::Section:: *)
(* error messages *)


getAE::incorrectNumberOfMomenta="The expression has more or less legs than provided by `1`. Adjust the number of legs.";

getAE::undefinedFields="Not all field(s) `1` are defined yet. Use defineFieldsSpecific to do so.";

getAE::incorrectFormat="The format of the list of external momenta and indices is not correct.
It has to be of the form {{field1, ind, mom, specific index 1, specific index 2, ...}, {field2, ...}, ...},
where field1 is a field, ind its index as given by the output of doDSE or doRGE, mom its momentum and specific index i are the indices of the field 1.
\nThe error was in `1`."; 

getAE::noFieldsGiven="No fields are indicated in the second argument of getAE, which has to be of the form {{field1, ind, mom, specific index 1, specific index 2, ...}, {field2, ...}, ...},
where field1 is a field, ind its index as given by the output of doDSE or doRGE, mom its momentum and specific index i are the indices of the field 1.
\nThe error was in `1`.";

getAE::fieldsAndIndicesDontMatch="The indices given for the external legs do not appear as indices of the indicated fields.
Make sure that the indices and fields of the external legs match.
The error was in `1`.";


(* :: Section:: *)
(* indices *)


putIndices[a_?NumericQ,___]:=a;

putIndices[a_Plus|a_Times,rest___]:=putIndices[#,rest]&/@a;

putIndices[a_List,rest___]/;Not@FreeQ[a,op,\[Infinity]]:=putIndices[#,rest]&/@a;

putIndices[a_op,indList_List,fields_List]:=Module[
	{fieldsInds},
fieldsInds=Cases[a,{_,_},2];
a/.List@@((#:>putIndices[#,indList,fields])&/@fieldsInds)
];

putIndices[{{Q_, q_},p_}, indList_List, fields_List]:= Module[{fermions, allIndList, relevantInds},
  
  (*fermions = Cases[fields, {_, _}];*)
  fermions=Select[fields,fermionQ];
  
  (* add the anti-fields to indList *)
  
  (*allIndList = 
   Join[indList, 
    Flatten[Cases[
        indList, {#[[1]], inds___} :> {#[[2]], inds}] & /@ 
      fermions, 1]];*)
  allIndList = Union[indList,indList/.{a_?fieldQ,inds___}:>{antiField[a],inds}];
   
  
  relevantInds = Rest@Flatten[Cases[allIndList, {Q, ___}]];
  (* put momentum p at first place in agreement with other functions, e.g. getFR *)
  Q[p,Sequence @@ (q[#] & /@ relevantInds)]
];


putIndices[{Q_, q_}, indList_List, fields_List]:=
 Module[{fermions, allIndList, relevantInds},
  
  (*fermions = Cases[fields, {_, _}];*)
  fermions=Select[fields,fermionQ];
  
  (* add the anti-fields to indList *)
  
 (* allIndList = 
   Join[indList, 
    Flatten[Cases[
        indList, {#[[1]], inds___} :> {#[[2]], inds}] & /@ 
      Cases[fields, {_, _}], 1]];*)
  allIndList = Union[indList,indList/.{a_?fieldQ,inds___}:>{antiField[a],inds}];
  relevantInds = Rest@Flatten[Cases[allIndList, {Q, ___}]];
  Q[Sequence @@ (q[#] & /@ relevantInds)]
  ];


addIndices[]:=Options[createDummyListUnique];
addIndices[{}]:=Options[createDummyListUnique];
(* several new indices *)
addIndices[{{name1_,dummies1_List},rest___}]:=Fold[addIndices[#2]&,{},{{name1,dummies1},rest}];
(* one new index *)
addIndices[{name_,dummies_List}]:=Module[
	{oldIndexDummies},
	
	(* known indices up to now *)
	oldIndexDummies=Options[createDummyListUnique];
	
	(* add new indices, overwrite old instances where necessary *)
	Options[createDummyListUnique] =  Join[{dummyNames[name] -> dummies},oldIndexDummies]
		/.{a___,dummyNames[b_]->c_List,d___,dummyNames[b_]->f_List,g___}:>{a,dummyNames[b]->c,d,g}
]


removeIndices[name_,rest__]:=removeIndices[{name,rest}];
removeIndices[names_List]:=Fold[removeIndices[#2]&,{},names];
removeIndices[name_]:=Options[createDummyListUnique] =  DeleteCases[Options[createDummyListUnique], dummyNames[name] -> _List]

resetIndices[]:=Options[createDummyListUnique]={dummyNames[adj]->{Global`a,Global`b,Global`c,Global`d,Global`e,Global`f,Global`g,Global`h,Global`i,Global`j},
	dummyNames[lor]->{Global`\[Mu],Global`\[Nu],Global`\[Rho],Global`\[Sigma],Global`\[Tau],Global`\[Alpha],Global`\[Beta],Global`\[Gamma],Global`\[Delta],Global`\[Epsilon]}};


(* create a unique list of dummies; note that in DoDSERGE another function createDummyList is used to alleviate the identification of diagrams *)
createDummyListUnique[a_Integer,type_,opts___]:=Module[{dNames},
	dNames=Join[{opts},type/.Options@createDummyListUnique];
	Take[Flatten@Table[Unique[dNames],{Ceiling[a/Length@dNames]}],a]
];




(* ::Section:: *)
(* momenta *)


(* add momentum to field and index *)
(* several momenta *)
addMomentum[b_, c_List, momenta_List] := 
  Fold[addMomentum[#1, #2[[1]], #2[[2]]] &, b, 
   Transpose[{c, momenta}]];

(* single momentum *)
addMomentum[b_, {Q_, q_}, mom_] := b /. {Q, q} :> {{Q, q}, mom};


(* momenta in props: for both indices the same and the attached \
vertices should also get the correct momentum *)

flowProps[b_, props_List] := 
  b /. Replace[
    List @@ # & /@ (Sort[#, (MatchQ[{#1, #2}, 
            List[{_?(Not@ListQ@# &), _}, {{_, _}, _}]] &)] & /@ 
       props), {{Q1_, q1_}, {{Q2_, q2_}, 
       p_}} :> ({Q1, q1} :> {{Q1, q1}, -p}), {1}];

insertMomenta[a_?NumericQ, rest___] := a;

insertMomenta[a_Times | a_Plus, rest___] := 
  insertMomenta[#, rest] & /@ a;

insertMomenta[a_List,rest___]:=insertMomenta[#,rest]&/@a;


(* no momenta to insert, i.e. an RGE vacuum diagram; only add internal momentum *)
insertMomenta[op[a___,P[b1_,b2_],c___],{}]/;Count[{a,c},dR[___],\[Infinity]]==1&&Count[{a,c},P[___],\[Infinity]]==0:=
    op[a,P[b1,b2],c] /. {b1 :> {b1, $loopMomenta[[1]]}, b2 :> {b2, -$loopMomenta[[1]]}};




(* treat tadpole extra *)

insertMomenta[op[S[a__], P[b1_, b2_]], extMomenta_List] := 
  Module[{external, totalOp},
   
   totalOp = op[S[a], P[b1, b2]];
   
   (* external indices *)
   
   external = 
    Cases[totalOp, {Q_, q_} /; 
      Count[totalOp, q, \[Infinity]] == 1, \[Infinity]];
   
   (* add external momenta and that for the loop *)

   addMomentum[
    totalOp /. {b1 :> {b1, $loopMomenta[[1]]}, 
      b2 :> {b2, -$loopMomenta[[1]]}},
    external, extMomenta]
   
   ];



(* treat expressions with only one entry in op, e.g., bare diagrams, extra *)

insertMomenta[op[SPV_[a__]], extMomenta_List]/;SPV==V||SPV==S||SPV==P||SPV==CO:= 
  Module[{external, totalOp},

   totalOp = op[SPV[a]];
   
   (* external indices *)
   
   external = 
    Cases[totalOp, {Q_, q_} /; 
      Count[totalOp, q, \[Infinity]] == 1, \[Infinity]];
   
   (* add external momenta and that for the loop *)
   
   addMomentum[totalOp, external, extMomenta]
   
   ];


insertMomenta[a_op, extMomenta_List] := 
 Module[{momentaToInsert, external, noListQ,  extAdded, props, regs,
   regRules, startVertex, startVertexIntLegs, startVertexExtLegs, loopMomenta, 
   newVertex, flowed,extFields,extFieldsAdded,extFieldsMomenta,
   verts, selfProps, momentaSelfProps, selfPropsAdded},
 
  momentaToInsert[b_, c_] := {c, -b - c};
  momentaToInsert[b_List, c_] := {Sequence @@ Most@b, Last@b - c, c};
  
  (* external indices *)
  
  external = 
   Cases[a, {Q_, q_} /; Count[a, q, \[Infinity]] == 1, \[Infinity]];
  
  noListQ[b_] := Not@ListQ@b;
  
  (* add momenta to external indices *)
  extAdded = addMomentum[a, external, extMomenta];
  
  extFields=Cases[a,{_?fieldQ,_}];

  (*extFieldsMomenta=Take[$loopMomenta,-Length@extFields];*)
  extFieldsMomenta=Table[0,{Length@extFields}];
  
  extFieldsAdded=addMomentum[extAdded, extFields, extFieldsMomenta];
  
  (* start at one vertex; take its indices *)
  
  startVertex = First@Cases[extFieldsAdded, S[b__] | V[b__] | CO[b__] :> {b}];
  (* the internal legs and the external one *)
  
  startVertexIntLegs = 
   Select[startVertex, FreeQ[{#}, Alternatives @@ Join[extMomenta,extFieldsMomenta]] &];
   (*Select[startVertex, FreeQ[{#}, Alternatives @@ extMomenta] &];*)
  
  startVertexExtLegs = Complement[startVertex, startVertexIntLegs];
  
  (* do self-loops *)
  (* get propagators and vertices and check which propagator arguments appear in the same vertex *)
  props = Cases[a, P[___]];
  (* get regulator insertions: to be treated special; actually not needed currently, but keep code for future use *)
  regs = Cases[a, dR[___]];
  (* rules to delete the dummy indices in P dR P *)
  regRules = (# -> Sequence[]) & /@ Flatten[List @@@ regs, 1];
  
  (* extract the propagators attached to the regulators *)
  regs = Cases[a, P[b___, {c_,Alternatives@@#[[All,2]]}, d___]]&/@regs;
  (* merge into one propagator representing the self-loop *)
  regs = regs /. regRules /. {P[b_List], P[c_List]} :> P[b, c];
  
  verts = Cases[a, V[___]|CO[___]|S[___]|dR[___]];
  (* get all propagators closing on the same vertex *)
  selfProps = Select[props, Position[verts,#[[1]]][[1, 1]]==Position[verts,#[[2]]][[1, 1]]&];
  
  (* and their momenta *)
  momentaSelfProps = $loopMomenta[[1;;Length@selfProps]];
  (* add those momenta *)
  selfPropsAdded = addMomentum[extFieldsAdded, Flatten[List@@@selfProps,1], Flatten[Transpose[{momentaSelfProps,-momentaSelfProps}],1]];
  
  (* remove the self-loops from the internal legs *)
  startVertexIntLegs = Complement[startVertexIntLegs, Flatten[List@@@selfProps,1]];
  
  (* determine the loop momenta *)
  
  If[Length[startVertexIntLegs]>1,
  	loopMomenta = 
   	Fold[momentaToInsert, (Plus @@ startVertexExtLegs)[[2]], 
   	$loopMomenta[[Length@selfProps+1;;Length@selfProps+Length@startVertexIntLegs-1]]],
   	(* no legs left, e.g., tadpoles *)
   	loopMomenta = {}
  ];
  
  (* replace the indices with the new momenta *)
  newVertex = addMomentum[selfPropsAdded, startVertexIntLegs, loopMomenta];
  
  (* the propagators with one new momentum *)
  props = Cases[newVertex, 
    P[{{_, _}, _}, {_?(noListQ@# &), _}] | 
     P[{_?(noListQ@# &), _}, {{_, _}, _}]];
  
  (* propagate the new momentum through the propagators and to the vertices on the other side *)
  flowed = flowProps[newVertex, props];

  FixedPoint[flowDiagram, flowed, 20]
  ];
  
  
(* nothing to do if everything is replaced *)

flowDiagram[a_op] /; 
   Cases[a, 
     P[{_?(Not@ListQ@# &), _}, {_?(Not@
           ListQ@# &), _}], \[Infinity]] == {} := a;
flowDiagram[a_op] := 
 Module[{startVertex, indWMomenta, indWOMomenta, newMomentum, newOp, 
   props, newLoopMomentum, usedLoopMomenta},
  
  (* choose one vertex that has only one leg without momentum *)
  (*startVertex = 
   First@Cases[a, 
     V[b__] | S[b__] | dR[b__] | CO[b__] :> {b} /; 
       Count[{b}, {_?(Not@ListQ@# &), _}] == 1];*)
  startVertex = 
   First@Cases[a,  
     V[b__] | S[b__] | dR[b__] | CO[b__] :> {b} /; 
       Count[{b}, {_?(Not@ListQ@# &), _}] >= 1];     
  (*indWMomenta = Cases[startVertex, {b_List, c_}];
  indWOMomenta = First@Complement[startVertex, indWMomenta];*)
  (* get legs with and without momenta *)
  indWMomenta = Cases[startVertex, {b_List, c_}];
  indWOMomenta = Complement[startVertex, indWMomenta];
  (*indWOMomenta = First@Cases[startVertex, {b_,c_}/;Head@b=!=List];
  indWMomenta = Complement[List@@startVertex, {indWOMomenta}];*)
  
    (* for exactly one leg without momentum continue by momentum conservation, otherwise add new loop momentum *)
  If[Length@indWOMomenta==1,
  	newLoopMomentum=-Plus @@ indWMomenta[[All, 2]],
  	
  	(* get an unused loop momentum *)
    usedLoopMomenta = Union@Cases[a, _?(MemberQ[$loopMomenta, #] &), Infinity];
    newLoopMomentum = First[Complement[$loopMomenta, usedLoopMomenta]];
    newLoopMomentum = First[$loopMomenta/.(#:>Sequence[]&/@usedLoopMomenta)];  	
  ];
  
  (* take first entry to put a momentum there *)
  indWOMomenta=First@indWOMomenta;
  
  (* determine the new momentum *)
  
  (*newMomentum = {indWOMomenta, -Plus @@ indWMomenta[[All, 2]]};*)
  newMomentum = {indWOMomenta, newLoopMomentum};
  
  (* add the new momentum *)
  
  newOp = addMomentum[a, newMomentum[[1]], newMomentum[[2]]];
  
  (* get the new propagator *)
  
  props = Cases[newOp, 
    P[{{_, _}, _}, {_?(Not@ListQ@# &), _}] | 
     P[{_?(Not@ListQ@# &), _}, {{_, _}, _}]];
  
  (* and propagate its momentum to the vertex at the other end *)
  
  flowProps[newOp, props]
  ];
  


(* ::Section:: *)
(* contract indices *)


(* contract indices using the function given, e.g. Contract or SUNSimplify of FeynCalc or a user function *)
(*contractIndices[a_?NumericQ,___]:=a;

contractIndices[a_Times|a_Plus,rest___]:=contractIndices[#,rest]&/@a;*)

contractIndices[a_,rest___]/;Not@FreeQ[a,op,\[Infinity]]:=contractIndices[a/.op[b__]:>Times[b],rest];

contractIndices[a_,contractFunctions_List,rest___]:=Fold[contractIndices[#1,#2,rest]&,a,contractFunctions];

contractIndices[a_,contractFunction_,opts___?OptionQ]:=contractFunction[a];



(* ::Sectionn:: *)
(* replace indices *)


(* replace the indices by "nice" expressions
replace all indices at once, i.e. do not resolve to the level of ops, so that external indices are the same everywhere *)


replaceIndices[a_,types_List,rest___]:=Fold[replaceIndices[#1,#2,rest]&,a,types];

replaceIndices[a_, type_, rest___]:=Module[
	{oldInds, indRules},
	
	(* get all indices in expression *)
	oldInds=Union@Cases[a, _[type], \[Infinity]];

	(* create rules to replace old by new indices *)
	indRules=Thread@Rule[oldInds, createDummyListUnique[Length@oldInds,dummyNames[type]]];
	
	a/.indRules
	
]

(* if a field has two indices of the same type the index will appear double in the field (four times in total for an internal index) -> rename the second one *)
removeDoubleIndices[exp_,fields_List]:=Module[{doubleIndexRules, allFields},
	
	allFields=Union[fields,antiField/@fields];

	(* rules that take rename the second index of double indices *)
	doubleIndexRules=Flatten[{op[
		h___, 
		vert_[d___, #[mom1_, otherInds1___, ind_[indtype_], otherInds2___, ind_[indtype_], otherInds3___], e___],
		b___,
		P[f___, #[mom2_, otherInd2___, ind_[indtype_], otherInds2___, ind_[indtype_], otherInds3___], g___],
		c___] :> Module[{dummy},
		
		(* get a appropriate name for the new index name *)
		dummy=createDummyListUnique[1,dummyNames[indtype]][[1]];

   		op[h, 
   			vert[d, #[mom1, otherInds1, ind[indtype], otherInds2, dummy[indtype], otherInds3], e],
   			b,
   			P[f, #[mom2, otherInds1, ind[indtype], otherInds2, dummy[indtype], otherInds3], g], c]
		],
		(*reverse order*)
		op[
		h___, 
		P[f___, #[mom2_, otherInd2___, ind_[indtype_], otherInds2___, ind_[indtype_], otherInds3___], g___],
		b___,
		vert_[d___, #[mom1_, otherInds1___, ind_[indtype_], otherInds2___, ind_[indtype_], otherInds3___], e___],
		c___] :> Module[{dummy},
		
		(* get a appropriate name for the new index name *)
		dummy=createDummyListUnique[1,dummyNames[indtype]][[1]];

   		op[h, 
   			vert[d, #[mom1, otherInds1, ind[indtype], otherInds2, dummy[indtype], otherInds3], e],
   			b,
   			P[f, #[mom2, otherInds1, ind[indtype], otherInds2, dummy[indtype], otherInds3], g], c]
		]}
     & /@ allFields
     ];
     
     exp//.doubleIndexRules
]



(* ::Section:: *)
(* convert symbolic expressions into algebraic ones *)


(* combine all functions necessary into one function; if demanded by the option save, save the result *)

(* if no contractFunctions are given (standard case), provide an empty list to use the original functions *)
getAE[a_, momenta_List,opts___?OptionQ]/;opts=!={}:=getAE[a, momenta, {}, opts];

(* fields are not defined for getAE *)
getAE[a_, momenta_List, contractFunctions_List, opts___?OptionQ]/;Not[
	And@@((ListQ@indices[#])&/@Flatten@momenta[[All,1]])]:=Message[getAE::undefinedFields,Union@momenta[[All,1]]];

(* in case the user has not defined an index yet, it is defined with standard index names, a,b,... *)
getAE[a_, momenta_List, contractFunctions_List, opts___?OptionQ]/;Not[
	And@@(MemberQ[Cases[Options[createDummyListUnique], dummyNames[b_]:>b, 2], #] & /@Union[Flatten[indices/@Flatten[Cases[a,_?fieldQ,\[Infinity]]]]])]:=Module[
	{undefinedInds, allInds,fields},
	fields=Cases[a,_?fieldQ,\[Infinity]];
	allInds=Union[Flatten[indices/@Flatten[fields]]];
	undefinedInds=Select[allInds,Not@MemberQ[Cases[Options[createDummyListUnique], dummyNames[c_]:>c, 2],#]&];
	addIndices[{#,{Global`a,Global`b,Global`c,Global`d,Global`e,Global`f,Global`g,Global`h}}]&/@undefinedInds;
	(* call the function itself again; the indices should now be defined *)
	getAE[a,momenta,contractFunctions,opts]
	]; 

(* in case there is not the right amount of momenta given *)
getAE[a_,  momenta_List, contractFunctions_List, opts___?OptionQ]/;
Count[Tally[Cases[Cases[Replace[a, op[b__] :> {op[b]}(*in case a=op[___]*)], op[___], \[Infinity]][[1]](*take first instance*), {_?fieldQ, _}, \[Infinity]]], {_, 1}]!=Length@momenta:=Message[getAE::incorrectNumberOfMomenta,momenta]; 

(* saving the results *)
getAE[a_,  extMomentaIndices_List, contractFunctions_List, opts___?OptionQ]/;(save/.Join[{opts},Options@getAE]):=
	getAE[a,  extMomentaIndices, contractFunctions, opts]=getAEDo[Expand@a, extMomentaIndices, contractFunctions, opts]

(* not the correct format of extMomentaIndices *)
getAE[a_,  extMomentaIndices_List, contractFunctions_List, opts___?OptionQ]/;
	Not[And@@(fieldQ/@extMomentaIndices[[All,1]])]:=Message[getAE::noFieldsGiven,extMomentaIndices]

(* external indices do not match the fields *)
getAE[a_,extMomentaIndices_List, rest___]/;Or@@(Cases[a,#,Infinity]==={}&/@extMomentaIndices[[All,{1,2}]]):=
	Message[getAE::fieldsAndIndicesDontMatch,extMomentaIndices];

(* add standard indices *)
getAE[a_,  extMomentaIndices_List, contractFunctions_List, opts___?OptionQ]/;
	Not[And@@(fieldQ/@extMomentaIndices[[All,1]])(*all first entries are fields*)]||
	Not[And@@(Length@#-3==Length@DoFun`DoFR`indices@#[[1]]&/@extMomentaIndices)](*numbers of given indices are correct*):=
	(Message[getAE::incorrectFormat,extMomentaIndices]
	(*getAE[a,Insert[#[[1]],#[[2]],2]&/@Transpose[{extMomentaIndices,Take[$externalIndices,Length@extMomentaIndices]}],contractFunctions,opts]*)
	);

getAE[a_,  extMomentaIndices_List, contractFunctions_List, opts___?OptionQ]:=
	getAEDo[Expand@a,  extMomentaIndices, contractFunctions, opts]

(* the vacuum diagrams of RGEs requires special treatment *)
(*getAEDo[a_,indicesList_List, fields_List, momenta_List, contractFunctions_List, opts___?OptionQ]:=*)


(* if only external momenta and no indices are given the momenta are brought into canonical form -> the program will leave the intermediary indices untouched;
mainly for compatibility reasons *)
getAEDo[a_, momenta_List, rest___]/;Depth[momenta]==2&&momenta=!={}:=
	getAEDo[a, {#}&/@momenta, rest];

getAEDo[a_, extMomentaIndices_List, contractFunctions_List, opts___?OptionQ]:=
	(*getAE[a,indicesList,  momenta, contractFunctions, opts]=*)Module[{fields,uniqueExt,legRules,indicesList,
	momentumRules,indsTypes,momentaInserted,indsPutIn,doubleIndsReplaced,indsReplaced,fullExpression, explicitTF,uniqueExtOrdered},
fields=Union@Cases[a,_?fieldQ,\[Infinity]];

(* get the fields and their indices *)
(*indicesList=fields2/. {Q1_?fermionQ, Q2_?antiFermionQ} :> Q1 /. Q3_?fieldQ :> {Q3, Sequence @@ indices@Q3};*)
indicesList=fields /. Q3_?fieldQ :> {Q3, Sequence @@ indices@Q3};

(* create unique external momenta which are replaced at the end *)
uniqueExt=Table[Unique["ext"],{Length@extMomentaIndices}];

(* rules for replacing indices and momenta in the external legs *)
(*legRules=Rule[Q_?fieldQ[#[[1]],___],Q[Sequence@@#[[2]]]]&/@Transpose[{uniqueExt,extMomentaIndices}];*)
legRules = Rule[#[[1]][#[[3]], #[[2]][___], ___], #[[1]][Sequence @@ Rest@Rest@#]] & /@ extMomentaIndices;

(* rules for replacing the unique external momenta by the user-given momenta *)
(*momentumRules=Thread[Rule[uniqueExt,extMomentaIndices[[All,1]]]];*)

(* should explicit expressions be used? *)
explicitTF=explicit/.Join[{opts},Options@getAE];

(* get the types of existing indices *)
indsTypes=Union@Flatten[Rest /@ indicesList];

(* transform to list: unique identification of graphs and expressions possible;
insert momenta *)
momentaInserted=insertMomenta[a/.Plus:>List, uniqueExt];
(* if there is only one graph, make it into a list *)
momentaInserted=Flatten@{momentaInserted};
(* derive for each graph the ordered list of unique external momenta corresponding to the order given in extMomentaIndices *) 
uniqueExtOrdered=Function[graph,(Flatten[Cases[graph, {{#[[1]], #[[2]]}, mom_}, \[Infinity]],1] & /@ extMomentaIndices)[[All,2]]]/@momentaInserted;

(* rules for replacing the unique external momenta by the user-given momenta for each single graph *)
momentumRules=Function[graph,Thread[Rule[graph,extMomentaIndices[[All,3]]]]]/@uniqueExtOrdered;

(* put in indices; complicated replacement of momenta (for each single graph) necessary to guarantee correspondence with given external momenta *)
(*indsPutIn=putIndices[momentaInserted, indicesList, fields]/.legRules/.momentumRules;*)
indsPutIn=MapThread[ReplaceAll,{putIndices[momentaInserted, indicesList, fields],momentumRules}]/.legRules;

doubleIndsReplaced=removeDoubleIndices[indsPutIn,fields];

(* replace indices by shorter expressions *)
indsReplaced=replaceIndices[doubleIndsReplaced, indsTypes];

(* transform to explicit expressions; therefor Feynman rules are needed in a special form;
remove op function *)
fullExpression=indsReplaced /. b_op :> Times @@ b /. {
  V[c___] :> V[c, explicit -> explicitTF], 
  CO[c___] :> CO[c, explicit -> explicitTF],
  S[c___] :> S[c, explicit -> explicitTF], 
  P[c___] :> P[c, explicit -> explicitTF],
  dR[c___] :> dR[c, explicit -> explicitTF]};

(* contract indices as wanted *)
contractIndices[fullExpression, contractFunctions]

]




(* ::Section:: *)
(* test functions for debugging *)


  
(* a test function for expressions already having momenta assigned (but the final expressions) *)
checkMomentumConservation[a_?NumericQ] := 1;
checkMomentumConservation[a_Times | a_Plus] := 
  checkMomentumConservation /@ a;
checkMomentumConservation[a_op] := 
 Apply[Plus, 
  Replace[Cases[a, V[b__] | S[b__] |CO[b__] :> {b}],{{Q_, q_}, p_} :> 
    p, {2}], {1}];


(* a test function for checking if any index appears more than twice *)

indicesTestAlg[a_]/;Not@FreeQ[a,{_,_}]:=Module[
{inds},

inds=Transpose[Cases[a,{_,ind_},Infinity]][[2]];
Union@Select[{#,Count[inds,#]}&/@inds,#[[2]]>2&][[All,1]]
];

indicesTestAlg[a_]:={};


checkIndicesAlg[a_]/;(b=Select[Cases[a,op[___],\[Infinity]],(indicesTestAlg[#]!={}&)])!={}:=
	Message[checkIndicesAlg::multipleIndices,b,indicesTestAlg[b]];

checkIndicesAlg[a_]:=Message[checkIndicesAlg::ok];




End[]

EndPackage[]
