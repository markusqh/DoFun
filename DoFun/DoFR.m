(* ::Package:: *)

(* Mathematica Package *)

(* published under GNU General Public License v3.0 *)




(* version history *)
(*  0.1 (Nov. 9, 2010): copying code from development notebook
    0.2: basic structure finished
    0.3 (Dec. 10, 2010): improved delta integrationd performance
    0.4 (Jan. 9, 2011)  (never released): consistent version for intensive field tests
    2.0.0 (25.2.2011): First public release of DoFun via http://www.tpi.uni-jena.de/qfphysics/homepage/mhub/DoFun/.
    2.0.1 (22.1.2013, never released): 
  		-) Added SyntaxInformation for fields in defineFieldsSpecific.
    2.0.2 (22.1.2015):
		-) Bugfix in integrateMomenta when momenta were only one character long.
		-) Bugfix in integrateMomenta when momentum is negative.
    2.0.3 (19.1.2017): no changes
    2.0.4 (5.12.2017): no changes
    3.0.0 (31.7.2019):
    	-) modified defineFieldsSpecific
    	-) modified derivF and getFR to left-derivatives
*)


BeginPackage["DoFun`DoFR`", { "DoFun`DoDSERGE`", "DoFun`DoAE`"}]


$DoFRVersion="3.0.0";

If[Not@FreeQ[Contexts[],"DoFun`"],DoFun`DoFR`$doFRStartMessage=False];
If[DoFun`DoFR`$doFRStartMessage=!=False,
	Print["\nPackage DoFR loaded.
\nVersion "<> $DoFRVersion <>
"\nJens Braun, Anton K. Cyrol, Markus Q. Huber, Jan. M. Pawlowski 2010-2019\n
\nDetails at https://github.com/markusqh/DoFun/."];
];



(* ::Section::Closed:: *)
(* usages *)

convertAction::usage="convertAction[ac] converts a given physical action ac into a form suitable for computation, i.e., with proper dummy indices and momenta.

Example:
setFields[{\[CurlyPhi]}];
defineFieldsSpecific[{\[CurlyPhi][momentum, type]}];
convertAction[ 1/2 (p^2 Z[p^2] + R[k]) op[\[CurlyPhi][p, i], \[CurlyPhi][-p,i]] + 
 \[Lambda]2/8 op[\[CurlyPhi][p, i], \[CurlyPhi][q, i], \[CurlyPhi][r, j], \[CurlyPhi][-p - q - r, j]]]
";

defineFieldsSpecific::usage="defineFieldsSpecific[fs] defines indices of fields fs.
The indices of a field can be obtained with indices[field].

Example:
setFields[{A,{c,cb}}];
defineFieldsSpecific[{A[momentum,adjoint,lorentz], c[momentum,adjoint], cb[momentum,adjoint]}]
indices/@{A,c,cb}
";

delta::usage="delta[a,b] represents an arbitrary Kronecker delta with indices a and b.
delta[ind,a,b] corresponds to a Kronecker delta for the index ind.

Examples:
{delta[1,1], delta[lorentz, mu, nu]}
";

deltam::usage="deltam[p1+p2] is the momentum delta distribution (2 \[Pi])^d \[Delta](p1+p2).
";

der::usage="der[expr, f] denotes the derivative of expr with respect to the field f.

Example:
der[U, phi[q, i]]
";

dim::usage="dim[ind] represents the dimension of the representation of the index ind. Values can be assigned by the user.

Example:
dim[adjoint]:=Nc;
integrateDeltas[delta[adjoint,a,b]delta[adjoint,b,a]]
";

getFR::usage="getFR[ac, fs] derives the Feynman rule for the n-point function with the legs fs from the action ac.

Example:
setFields[{\[CurlyPhi]}]
defineFieldsSpecific[{\[CurlyPhi][mom, type]}];
getFR[convertAction[1/2 p2 op[\[CurlyPhi][q1, j], \[CurlyPhi][-q1, j]] + 
 1/8 \[Lambda] op[\[CurlyPhi][q1, j], \[CurlyPhi][q2, j], \[CurlyPhi][q3, l], \[CurlyPhi][-q1 - q2 - q3, l]]], 
{\[CurlyPhi][p1,i],\[CurlyPhi][p2,j]}]
";

indices::usage="indices[f] gives the types of indices of a field f.

Example:
setFields[{phi}];
defineFieldsSpecific[{phi[momentum,indexOfPhi]}];
indices[phi]
";

integrateDeltas::usage="integrateDeltas[expr] contracts indices of Kronecker deltas in expr.

Examples:
integrateDeltas[delta[a, b] delta[b, c]]
integrateDeltas[delta[ind1, a, b] delta[ind1, a, b]]
integrateDeltas[delta[ind1, a, b] delta[ind1, b, c]]
";

integrateMomenta::usage="integrateMomenta[expr] integrates out internal momenta in expr, denoted by q$i, where i is a running number, in momentum delta distributions.
integrateMomenta[expr, mom] integrates out the momenta mom in expr. mom can be a single momentum or a list of momenta.

Examples:
integrateMomenta[deltam[p1+q$1]deltam[q$1-p3]]
integrateMomenta[deltam[r1+p2]deltam[p2-r3],p2]
";

U::usage="U[f] represents a potential depending on f.
";




(* ::Section::Closed:: *)
(* options *)
(* define all functions which can contain indices *)


Options[integrateDeltas]:={delta}





Begin["`Private`"];

(* private functions (alphabetic)
	diffFields
	getDeltasForInt
	getDeltasForMixed
	getIndices
	integrateDeltaDo
	
*)


(* ::Section::Closed:: *)
(* Messages *)


convertAction::undefinedField="There is an undefined field or a field without specified indices in the expression. Use defineFieldsSpecific to define this field properly.
\nThe field(s) in this expression is/are `1`."

getFR::undefinedField=convertAction::undefinedField;

integrateDeltas::tooManyIndicesInDeltas="There is at least one index that appears more than twice in the deltas.
\nThe error is in the expression `1`."

integrateDeltas::sameIndexForDifferentTypes="There is at least once the same index for two different types of indices, e.g., delta[ind1,a,b]delta[ind2,b,c]. This is not supported by DoFR.
\nThe error is in expression `1`."

defineFieldsSpecific::fieldsNotSet="Not all fields in `1` are defined. Do so with setFields."




(* ::Section::Closed:: *)
(* miscellaneous functions *)


(* derivative of generic expressions exp; e.g. exp could be a potential U *)
der[der[exp_,field1_],field2_]:=der[exp,Flatten@{field1,field2}]


(* defining the fields *)

(* check if fields are defined *)
defineFieldsSpecific[fields_List, opts___?OptionQ] /; Not[And@@(fieldQ/@Flatten[fields][[All,0]])] := Message[defineFieldsSpecific::fieldsNotSet, fields];

defineFieldsSpecific[fields_List, opts___?OptionQ] := 
  Module[{},
  
   (* set upvalues: indices of fields *)   
   (Evaluate[#[[0]]] /: indices[Evaluate[#[[0]]]] := 
       List @@ Rest@#) & /@ Flatten@fields;
       
   (* define syntax information such that the user recognizes errors in input *)
   (SyntaxInformation[#[[0]]] = {"ArgumentsPattern" -> 
      Table[_, {Length@#}]}) & /@ Flatten@fields;
   
   fields
];
   

(* convert the "naive" user action into an expression with real dummy variables *)
(* error if an op expression contains undefined fields or fields with no indices specified *)
convertAction[a__ exp_op|exp_op,___]/;Module[{fields},fields=List @@ exp[[All,0]];Not[And@@(fieldQ /@ fields)]||And@@(indices[#] === # & /@ fields)]  :=(Message[convertAction::undefinedField,List @@ exp[[All,0]]];Abort[]);
convertAction[exp_,___]/;Module[{fields},fields=Union@Flatten[List @@@ Cases[exp, op[__], \[Infinity]][[All, All, 0]]];
	fields=!={}&&(Not[And@@(fieldQ /@ fields)]||And@@(indices[#] === # & /@ fields))] :=(Message[convertAction::undefinedField,Union[Flatten[List @@@ Cases[exp, op[__], \[Infinity]][[All, All, 0]]]]]);

convertAction[ac_Plus] := (convertAction[#] & /@ ac)

convertAction[ac_] /;Head@Expand@ac==Plus:= (convertAction[#] & /@ Expand@ac)

convertAction[ac_] := 
 Module[{fields,inds, newIndices, loopMoms, newLoopMoms, acExp},
  acExp=Expand@ac;
  fields=Cases[acExp,_?(fieldQ@#[[0]]&),\[Infinity]];
  inds = Union[Flatten[List @@ (List @@@ Rest /@ fields)]];
  newIndices = Table[insDummy[], {Length@inds}];
  loopMoms = 
   Union[Flatten[(First /@ List @@ (fields /. Plus :> List /. -q_ :> q))]];
  newLoopMoms = Table[Unique[Global`q],{Length@loopMoms}](*Take[DoFun`DoAE`$loopMomenta, Length@loopMoms]*);

  acExp /. Thread[Rule[inds, newIndices]] /. 
   Thread[Rule[loopMoms, newLoopMoms]]
  ]
  



(* ::Section::Closed:: *)
(* differentiation functions *)


derivF[lag_Plus, field_] := (derivF[#, field] & /@ lag)

derivF[lag_Times, field_] := (lag/# derivF[#, field]) & /@ Plus@@lag

derivF[exp_op, field_]/;FreeQ[exp,Evaluate[field[[0]]]]:=0
derivF[exp_op, field_] := Module[{fieldsInExp,orderedExps,moveField},

	(* moveField moves Grassmann fields acted upon by a derivative to the utmost left;
	it is always assumed that derivatives with respect to Grassmann fields act from the left *)	
	moveField[lExp_,lField_]/;grassmannQ@Head@lField:={lExp//.{op[a___,lField,b_?(grassmannQ@Head@# &),c___]:>-op[a,b,lField,c],
		op[a___,lField,b_?(cFieldQ@Head@# &),c___]:>op[a,b,lField,c]},lField};
	moveField[lExp_,lField_]:={lExp,lField};		
	
  fieldsInExp = List @@ Select[exp, Head@# == field[[0]] &];
  
  orderedExps=moveField[exp,#]&/@fieldsInExp;

  Plus@@Function[single,(deltam[single[[2,1]], field[[1]]] Times@@Thread[delta[indices[single[[2,0]]],List@@Rest@single[[2]], List@@Rest@field]] single[[1]]  /. single[[2]] :> Sequence[] )]/@orderedExps

  ]


derivF[exp_,field_]/;Not@FreeQ[exp,Head@field]:=der[exp,field]

derivF[exp_, field_] := 0

diffFields[exp_,fields_List]:=Fold[derivF[Expand@#1,#2]&,exp,fields]/. op[] :> 1




(* ::Section::Closed:: *)
(* momentum deltas *)
(* canonical ordering of momenta in delta distribution *)


deltam[p_, q_] := deltam[p - q]
deltam[-p_ - q_] := deltam[p + q]
deltam[q_ - p_] /; Sort[{p, q}] =!= {p, q} := deltam[-q + p]

(* public function for integrating; note that no internal momenta are integrated out, if a momentum for integration is given by the user *)
integrateMomenta[exp_Plus,rest___] := integrateMomenta[#,rest]& /@ exp

integrateMomenta[exp_Times,rest___]:=FixedPoint[integrateMomentum[#1,rest]&,exp,20]

integrateMomenta[exp_]:=exp

(* if no momenta are given explicitly, the internal momenta q$i are integrated out *)
integrateMomentum[exp_Times]/;Count[exp,deltam[___],\[Infinity]]>=1(*||(Count[exp,deltam[___],\[Infinity]]==1&&Count[exp,op[___],\[Infinity]]==1)*):=Module[
	{intMomenta,allDeltas},
	allDeltas=Cases[exp,deltam[___],\[Infinity]];
	intMomenta=Union[ToExpression/@Select[SymbolName /@ Cases[allDeltas, s_Symbol|-s_Symbol:>s, {3}]/.a_String:>Sequence[]/;StringLength@a==1(*drop too short expressions*), StringTake[#, 2] == "q$" &]];

	integrateMomentum[exp,intMomenta]	
]

integrateMomentum[exp_Times,momenta_List]:=Fold[integrateMomentum[#1,#2]&,exp,momenta];

integrateMomentum[exp_Times,momentum_]/;FreeQ[Cases[exp,deltam[___],\[Infinity]],momentum]:=exp;

(* integration routine where momentum is integrated out *)
integrateMomentum[exp_Times|exp_Power,momentum_]:=Module[
	{deltaToIntegrate,signOfMom=1},
	
	(* which delta to use for the integration *)
	deltaToIntegrate=First@Cases[exp, deltam[a__]/;Not@FreeQ[a,momentum], \[Infinity]];
	
	(* determine sign in front of momentum *)
	If[Not@FreeQ[deltaToIntegrate, -momentum], signOfMom = -1];

	(exp/deltaToIntegrate)/.momentum:>-signOfMom(Plus@@deltaToIntegrate)+momentum
]

integrateMomentum[exp_,rest___]:=exp



(* ::Section::Closed:: *)
(* generic deltas and their contractions *)


(* trace of a delta *)
delta[ind_,a_Symbol|a_dummy,a_Symbol|a_dummy]=dim[ind];
(* numeric values for delta *)
delta[ind___,a_Integer,a_Integer]:=1;
delta[ind___,a_Integer,b_Integer]/;b=!=c:=0;


(* do contraction of indices in delta *)
(* no integration to do; end of iteration *)
integrateDeltaDo[
  exp_, {{}}|{}] := exp (* to stop the fixed point iteration *)

(* for fixed point interations lists are used *)
integrateDeltaDo[exp_, de_List] := (exp/de[[1]]) /. Rule @@ de[[1]]

(* integrate delta *)
integrateDeltaDo[exp_, de_] := (exp/de) /. Rule @@ Take[de,-2];


(* get indices in an expression *)
getIndices[exp_] := 
 Cases[exp, 
  Alternatives @@ (#[___] & /@ 
     Options@integrateDeltas), \[Infinity]]


(* get mixed and internal deltas *)

getDeltasForMixed[exp_] := 
 Module[{expMod, mixedDeltas, intDeltaInds, allDeltas, extDeltaInds,allDeltaInds},
  (* disintuish between external and internal indices -> 
  three different deltas: two internal, 
  two external and mixed indices *)
  
  (* give every index its name so there are no mix-ups possible between the same name but two different indices *)
  expMod=exp/.delta[ind_,a_,b_]:>delta[{ind,a},{ind,b}];
  
  (* only take non-numeric indices; numeric indices shoul be treated with integrateNumDeltas *)
  allDeltas = Cases[exp, delta[___,_Symbol|_dummy,_Symbol|_dummy], \[Infinity]];
  allDeltaInds = Join[Sequence @@@ (allDeltas/.delta[ind_,a_,b_]:>delta[a,b])];
  extDeltaInds = 
   Select[allDeltaInds, Count[exp, #, \[Infinity]] == 1 &];
    
  intDeltaInds = 
   Union@Select[allDeltaInds, Count[exp, #, \[Infinity]] == 2 &];

   mixedDeltas = 
   Select[allDeltas, 
     MemberQ[intDeltaInds, Alternatives @@ (#/.delta[ind_,a_,b_]:>delta[a,b])] && 
       MemberQ[extDeltaInds, Alternatives @@ (#/.delta[ind_,a_,b_]:>delta[a,b])] &];
  
  {mixedDeltas, intDeltaInds}
]


getDeltasForInt[exp_] := 
 Module[{intDeltas,allDeltas,intInds,allInds(*allDeltas, intDeltas, extDeltas, allDeltaInds, intDeltaInds, 
   extDeltaInds, mixedDeltas, allInds, intInds*)},
  (* disintuish between external and internal indices -> 
  three different deltas: two internal, 
  two external and mixed indices *)
  
  (* dummy[...] is not a symbol, so take it into account explicitly *)
  allDeltas = Cases[exp, delta[___,_Symbol|_dummy,_Symbol|_dummy], \[Infinity]];
  allInds = Flatten[List @@@ getIndices[exp/.delta[ind_,a_,b_]:>delta[a,b]]](*Join[Sequence @@@ allDeltas]*);
  (*allInds = Flatten[List @@@ getIndices[exp]](*Join[Sequence @@@ allDeltas]*);*)

  intInds = 
   Union@Select[allInds, Count[exp, #, \[Infinity]] == 2 &];
  
  (* the {} at the end is to terminate the series of replacements; count from the end to avoid the index type *)
  intDeltas = 
   Append[Select[allDeltas, 
    MemberQ[intInds, #[[-1]]] && 
      MemberQ[intInds, #[[-2]]] &&#[[-1]]=!=#[[-2]] &],{}]
]


(* deltas with numeric values *)
integrateNumDeltas[exp_Plus] := integrateNumDeltas[#] & /@ exp;

integrateNumDeltas[exp_]:=Module[
	{allDeltas, allDeltaInds,intDeltaInds,allNumDeltasInt,expMod,squareRules},
	
 	(* trivial numerical delta replacements *)
 	squareRules={delta[ind___,a_Integer,b_dummy|b_Symbol]^2:>1,delta[ind___,b_dummy|b_Symbol,a_Integer]^2:>1};
 	
 	FixedPoint[Function[x,
 	expMod=x/.squareRules;
 	
	(* all deltas and their non-numeric indices *)
	allDeltas = Cases[expMod, delta[___], \[Infinity]];
 	allDeltaInds = Cases[Join[Sequence @@@ allDeltas],_?(Not@IntegerQ@#&)];
 	
 	(* determine all internal indices *)
 	intDeltaInds = Union@Select[Union@allDeltaInds, Count[x, #, \[Infinity]] >= 2 (*>= if the index appears also elsewhere, e.g., in the numerator of the propagator*)&];
	
 	(* get all deltas with one numeric value and one internal index *)
 	allNumDeltasInt=Cases[allDeltas,delta[ind___,a_?(MemberQ[intDeltaInds,#]&),b_Integer]|delta[ind___, b_Integer,a_?(MemberQ[intDeltaInds,#]&)]:>delta[a,b]];
 	
 	Fold[integrateDeltaDo[#1,#2]&,x,allNumDeltasInt]/.squareRules],Expand@exp,20]
 	
 	
]


(* contract the indices in deltas *)
integrateDeltas[exp_Plus|exp_List] := integrateDeltas[#] & /@ exp
integrateDeltas[exp_] /;(Head@Expand@exp)==Plus:=integrateDeltas[#] & /@ Expand@exp

(* error handler if any non-numeric index appears more often than twice in the deltas *)
integrateDeltas[exp_]/;Module[{allDeltaInds},
		(* all non-numeric indices; structure delta[d___,e_,c_] is for discarding the index in the delta *)
		allDeltaInds=Select[Cases[List@@exp/.delta[d___,e_,c_]^n_Integer:>Table[delta[e,c],{n}](*deal with powers->make lists*), delta[a__,f_,g_]:>Sequence[f,g], \[Infinity]],Not[NumericQ@#]&];

		(* check if any comes more than twice *)		
		Max[Count[allDeltaInds,#]&/@allDeltaInds]>2
]:=(Message[integrateDeltas::tooManyIndicesInDeltas,exp];Abort[];)

(* error handler if any non-numeric index appears in deltas for different index types *)
integrateDeltas[exp_]/;Module[{allDeltaInds},
		(* all non-numeric indices; gives lists of all indices together with their index type *)
		allDeltaInds=Select[Cases[List@@exp, delta[a__,f_Symbol|f_dummy,g_Symbol|g_dummy]:>Sequence[{a,f},{a,g}], \[Infinity]],Not[NumericQ@#]&];
		
		(* gather indices and see if any appears for more than one index type *)		
		Max[Length /@ Union /@ GatherBy[allDeltaInds, #[[2]] &]]>1
]:=(Message[integrateDeltas::sameIndexForDifferentTypes,exp];Abort[];)


integrateDeltas[exp_] /;(Head@Expand@exp)==Times||(Head@Expand@exp)==Power:= 
 Module[{intDeltas, mixedDeltas,    intDeltaInds, modExp, mixedDeltasSorted,intIntegrated,mixedDeltasUse,intDeltasUse,squareRules},
      
     (* treat squares of deltas special *)
     squareRules={delta[ind___,a_Symbol|a_dummy,b_Symbol|b_dummy]^2:>delta[ind,a,a]};
     (* rules for delta[a,a] *)
     
     (* integrate out the internal deltas *)
     intIntegrated=FixedPoint[
     (
     intDeltas = getDeltasForInt[#];
     (* only use deltas that do not interfer with each other, i.e. they have no indices in common that are transformed simultaneously *)
     intDeltasUse=Union[Most@intDeltas, SameTest -> (Length@Union[{Sequence @@ Take[#1,-2], Sequence @@ Take[#2,-2]}] != 4 (*&& Length@#1 == 2 && Length@#2 == 2*) &)];
    
     (* sort the mixed Deltas in a canonical order: int, ext *)
     Fold[integrateDeltaDo[#1,#2]/.squareRules&,#(*from FixedPoint*),intDeltasUse]) &,integrateNumDeltas@exp, 20];
          
     {mixedDeltas, intDeltaInds} = getDeltasForMixed[intIntegrated];

     (* sort the mixed Deltas in a canonical order so they can be transformed into rules: delta[int, ext] *)
     mixedDeltasSorted=mixedDeltas/.delta[ind___,a_,b_]:>Join[delta[ind],Sort[delta[a,b], MemberQ[intDeltaInds, #] &]] ; 
     modExp = intIntegrated /. Thread@Rule[mixedDeltas, mixedDeltasSorted];
     
     (* take each internal index only once; count from the back to avoid possible index types at position 1 *) 
     mixedDeltasUse=Union[mixedDeltasSorted, SameTest -> (#1[[-2]] === #2[[-2]] &)];
          
     (* integrate all mixed Deltas and do some trivial numerical deltas *) 
     integrateNumDeltas@Fold[integrateDeltaDo[#1,#2]/.squareRules&,modExp,mixedDeltasUse]/.squareRules
  ]

integrateDeltas[exp_]:=exp



(* ::Section::Closed:: *)
(* function to derive the Feynman rules *)


(* fields not defined *)
getFR[action_,fields_List,rest___]/;Not[And@@(fieldQ /@ Union@fields[[All,0]])]||And@@(indices[#] === # & /@ Union@fields[[All,0]]):=(Message[getFR::undefinedField,Union@fields[[All,0]]]);

(* several FM rules *)
getFR[action_,{fields1_List,fields2___},opts___?OptionQ]:=getFR[action,#,opts]&/@{fields1,fields2}

getFR[action_,fields_List,opts___?OptionQ]:=Module[
	{nPoint, filtered, sign, orderedFields, counter=0},
	
	(* for vertices add a minus sign according to convention *)
	sign=Which[Length@fields>2,DoFun`DoDSERGE`$signConvention (-1),True,(+1)];

	
	(* filter out the terms with the appropriate fields to avoid unnecessary computations; keep everything in the broken phase *)
	filtered=If[(symmetry/.Join[{opts},Options@doRGE])==="broken"||Head@action==Plus (* in case there is only one term *),action,Plus@@Cases[action,(a___ op[b___]/;Sort[{b}[[All,0]]]==Sort[fields[[All,0]]])| U[__] | a__ U[__]],action];

    (* order anti-fermions due to the convention that fermion derivatives act from the right, i.e.,
       left before right fermions, and anti-fermions act from the left, i.e., right before left anti-fermions *)
    (*orderedFields=fields//.{{a___,b_?(fermionQ@Head@#&), c_?(antiFermionQ@Head@#&),d___}:>(counter++;{a,c,b,d}),
    	{a___,b_?(fermionQ@Head@#&), c_?(bosonQ@Head@#&),d___}:>{a,c,b,d},{a___,b_?(bosonQ@Head@#&), c_?(antiFermionQ@Head@#&),d___}:>{a,c,b,d}};*)
    (* DoFun 3.0: no order changed here, only left derivatives --> put Grassmann fields of derivatives in inverse order than in vertex,
    e.g., ghost-gluon vertex A cb c: fields=A, c, cb *)
    orderedFields=fields;
    
    nPoint=diffFields[filtered, orderedFields];

 	(* set fields to physical value, i.e.,to zero in the symmetric phase and take the remaining fields as vacuum fields in the broken phase *)
 	sign (-1)^counter If[(symmetry/.Join[{opts},Options@doRGE])==="broken",
 			integrateMomenta@integrateDeltas[nPoint]/.op[a__]:>Times@@{a},
 			integrateMomenta@integrateDeltas[nPoint]/.op[___]:>0,
 			integrateMomenta@integrateDeltas[nPoint]/.op[a__]:>Times@@{a}
 		
 	]
]



End[];

EndPackage[]
