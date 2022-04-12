(* ::Package:: *)

(* Mathematica Package *)

(* published under GNU General Public License v3.0 *)




(* ::Section:: *)
(* History *)


(* version history *)
(*  0.1: first running version
    0.2: performance improvements for derivation of higher vertex functions
   		no longer rules for the propagator necessary from the user, but still possible
    0.3: signs for fermion loops
    0.4: further simplifications in user input, list of interactions instead of action now sufficient
   
    1 (15.8.08): 1st published version on arXiv
    1.1 (16.09.08): -) standard vertex test added for fields with A <-> -A symmetry in action
        -) user does no longer have to specify indices in doDSE
        -) setSourcesZero: new argument interactions
        -) generateAction: user can specify bosonic and fermionic fields manually
        -) new options of doDSE: specificFieldDefinitions (used for mixed bosonic propagators in the action)
    1.2 (13.10.08): -) fixed bug when using doGrassmannTest->False and a user vertex test function
        -) new options indexStyle and factorStyle for DSEPlot
        -) (old) DSEPlot renamed to DSEPlotList;
           introduced DSEPlotGrid;
           (new) DSEPlot draws now complete DSEs including the left-hand side and prefactors
    1.2.1 (22.10.08): -) Exporting DSEPlots works now fine (the "" are not shown any longer)
    1.2.2 (27.10.08): -) added -1 exponents for propagators in DSEPlots
    1.2.3 (09.02.10): -) fixed bug: additional check for allowed propagators in setSourcesZero to get rid of wrong tadpole diagrams (they have no dressed vertex and thereby could circumvent the standard test
   		-) introduced variable $doDSEVersion
   		possibly there is a problem if you define cb,c or c,cb in the interaction list for the bare c propagator
   		-) introduced DoDSE`$startMessage to allow suppression of the message when loading the package
		-) added a property to op that is similar to the attribute Flat
		-) fixed bug for mixed complex propagators in getFermionList
		-) added Expand at several points to avoid problems with parentheses
		-) removed a correction factor 0.5 for the font size of the exponent on the lhs
    1.2.4 (never released): -) removed further correction factors 0.5 in fontsize of exponents (size looks ok now in pdfs compared to other font sizes)
   	    -) new functions: RGEPlot, doRGE
   	    -) right now DoDSE should not be loaded when loading DoRGE because the same context is used and wrong defitions will be used;
   	       this will be fixed when merging DoDSE and DoRGE
   	    -) new functions regulatorBox, regulatorCross
   	    -) new option regulatorSymbol
   	    -) new function grassmannQ
   	    -) new internal function changeOrder
   	    -) new function derivRGE
   	    -) new function getLoopNumber
   	    -) plotting of an equation even for a single graph is possible
   	    -) bugfix: indexStyle works for external fields
   	    -) new functions replaceVerticesBrokenPhase, newVertex, newVertexSum, update of setSourcesZero: allow the treatment of a broken phase, where the vacuum expectation values of fields do not vanish
   	    -) new option arrowHeadSize for DSEPlot
   	    -) new function derivPropagatorRGE
   	    -) fixed a bug in generateAction for theory with only one field which is fermionic
   	    -) new possibility to enter action via multiplicity of fields
   	    -) bugfix: specificFieldDefinitions did not work
    2.0.0 (25.2.2011): First public release of DoFun via http://www.tpi.uni-jena.de/qfphysics/homepage/mhub/DoFun/.
    2.0.1 (22.1.2013, never released):
		-) bugfix: in DSEPlotList with plotRules: Misidentification of fields like phiR as phi was possible. As a result the colors were not correct in the plots.
		-) new: options in DSEPlot and RGEPlot are forwarded to Grid
		-) bugfix in DSEPlot/RGEPlot: added Expand in condition to override possible simplifications
		-) DSEPlot/RGEPlot: default ImageSize 100 to avoid spacing problems with several rows
    2.0.2 (22.1.2015): no changes
    2.0.3 (19.1.2017):
    	-) fixed bug in doRGE: identifyGraphsRGE identified diagrams with opposite sign but different ordering of inner fermion legs. Introduced a canonical order for inner fermion legs to avoid this problem.
    2.0.4 (5.12.2017):
    	-) fixed bug in doRGE: Bugfix from 2.0.3 extended such that external fields are not affected thereby fixing a similar bug where a diagram was lost by identification due to sign problems.
    	-) Introduced a canonical ordering within fermion and antifermion fields.
    3.0.0 (31.7.2019):
    	-) generalization of diagram identifucation
    	-) basic rewriting of derivations of DSEs and RGEs to take into account Grassmann fields properly
	-) left derivatives instead of left/right derivatives
    	-) modifications of plotting (using Graph now)
    	-) derivation of equations for correlation functions of composite operators
        -) new functions:
        	* sortCanonical
        	* doCO
        	* COPlot
        	* getSigns
        	* sf
        	* extractDiagramType
        	* getDiagramType
        	* getVertexNumbers
        	* groupDiagrams
        	* cFieldQ
        	* disconnectedQ
        	* connectedQ
        	* getDisconnected
        	* getConnected
        	* onePIQ
        	* get1PI
        	* getNon1PI
        -) fields have a type now which is set by setFields and not guessed from input
        -) modifications in plotting: standalone propagators can be plotted now
        -) declared orderFermions deprecated; superseded by sortCanonical
        -) replaced identifyGraphs by identifyGraphsRGE --> only one function now which is called identifyGraphs; removed compareGraphs and compareGraphs2
        -) new options of DSEPlot/RGEPlot: bareVertexSymbol, vertexSymbol, coSymbol
        -) new symbols for plotting: diskSymbol, diskOpenSymbol, diskTinySymbol, triangleSymbol
        -) removed getExtGrassmannOrder and related functions
        -) removed specificFieldDefinitions
	3.0.1:
		-) introduced warning in sortCanonical when dummy fields are contained and the expression is not sorted
		-) introduced error handler when user tries to use the dummy field as field in setFields
	develop:
		-) grassmannTest is no longer done for RGEs, since ansatz for action fully determines the allowed vertices.
			Changing this behavior avoids problems with mixed fermionic propagators for RGEs. 
		-) replacedField for DSEs now handles sf expressions.
		-) $dummyField is no longer grassmann false but undefined to avoid wrong signs from sf.
*)



(* ::Section:: *)
(* Usages *)


BeginPackage["DoFun`DoDSERGE`"]

(* Exported symbols added here with SymbolName::usage *)

(* version of DoDSERGE *)
$doDSERGEVersion="3.0.1";

If[Not@FreeQ[Contexts[],"DoFun`"],DoFun`DoDSERGE`$doDSERGEStartMessage=False];
If[DoFun`DoDSERGE`$doDSERGEStartMessage=!=False,
	Print["Package DoDSERGE loaded.
\nVersion "<> $doDSERGEVersion <>
"\nJens Braun, Anton K. Cyrol, Reinhard Alkofer, Markus Q. Huber, Jan. M. Pawlowski, Kai Schwenzer, 2008-2019\n
\nDetails at https://github.com/markusqh/DoFun/."];
];




(* ::Section:: *)
(* Usages *)


(* symbols *)
 	
$bareVertexSymbol::usage="$bareVertexSymbol is the symbol for a bare vertex when using shortExpression.
Default value: S.
";

$compOpSymbol::usage="$compOpSymbol is the symbol for a composite operator when using shortExpression.
Default value: \[Omicron].";

$diagramTypes::usage="$diagramTypes is a list of all known diagram types.
";

$dummyField::usage="$dummyField is the super field representing all possible fields.
Default value: \[Phi].
";

$externalIndices::usage="$externalIndices contains the default names of external indices when none are supplied to doCO, doDSE or doRGE.
";

$fieldTypes::usage="$fieldTypes is the list of all existing field types.
";

$propagatorSymbol::usage="$propagatorSymbol is the symbol for a propagator when using shortExpression.
Default value: \[CapitalDelta].";

$regulatorInsertionSymbol::usage="$regulatorInsertionSymbol is the symbol for a regulator insertion when using shortExpression.
Default value: R.
";

$signConvention::usage="$signConvention sets the convention for the definition of vertices.
$signConvention=-1 (default since DoFun 3.0.0) means that vertices are the positive derivative of the effective action at the physical values of the fields.
$signConvention=1 (default before DoFun 3.0.0) means that vertices are the negative derivative of the effective action at the physical values of the fields.
";

$vertexSymbol::usage="$vertexSymbol is the symbol representing a vertex when using shortExpression.
Default value: \[CapitalGamma].
";

ansatz::usage="ansatz is an option of doDSE. It specifies which vertices are allowed in form of a list of possible interactions.
Note: For doRGE, the action corresponds to the ansatz for the effective average action.
See ?generateAction for possibilities on specifying interactions.
";

antiComplex::usage="antiComplex is the field type of a complex bosonic anti-field.
Properties of fields need to be set by setFields.
";

antiComplexFieldQ::usage="antiComplexFieldQ[expr] gives True if expr is a complex bosonic anti-field.  
";

antiFermion::usage="antiFermion is the field type of a Grassmann anti-field.
Properties of fields need to be set by setFields.
";

antiFermionQ::usage="antiFermionQ[expr] gives True if expr is a Grassmann anti-field.
";

antiField::usage="antiField[field] gives the anti-field of field.
Example:
setFields[{A}, {{c,cb}}, {{psi,psib}}];
antiField/@{A, c, cb, psi, psib}
";

bareVertexSymbol::usage="bareVertexSymbol is an option of COPlot, DSEPlot, DSEPlotList and RGEPlot. It determines how to draw bare vertices.

Possible values: boxSymbol, diskSymbol, triangleSymbol, diskTinySymbol, diskOpenSymbol, crossSymbol or a user-defined function which takes the coordinate of the regulator insertion as input.
Default value: diskTinySymbol.

Example:
setFields[{phi}, {}, {}]; DSEPlot[
 1/2 op[P[{phi, t1}, {phi, v1}], 
   S[{phi, i}, {phi, j}, {phi, v1}, {phi, t1}]], {{phi, Black}}, 
 bareVertexSymbol -> triangleSymbol]
";

boson::usage="boson is the field type of a real bosonic field.
Properties of fields need to be set by setFields.
";

bosonQ::usage="bosonQ[expr] gives True if expr is real bosonic field.
";

boxSymbol::usage="boxSymbol is a box graphic used for bareVertexSymbol, coSymbol, regulatorSymbol or vertexSymbol.
";

broken::usage="broken specifies for doDSE and doRGE whether a symmetry is intact or broken. In the latter case, the physical fields do not vanish.
";

cFieldQ::usage="cFieldQ[expr] gives True if expr is a commuting field.
";

checkAction::usage="checkAction[action] checks indices in action, i.e., it looks for free indices which should be absent for a properly defined action. Performs also checkSyntax and checkFields.\n

Example:
checkAction[op[S[{A, i1}, {B, i2}], {A, i1}, {B, j1}]]
";

checkAll::usage="checkAll[expr] performs a series of checks (checkIndices, checkSyntax, checkFields) on expr.\n

Example:
checkAll[op[ S[{A, i1}, {A, j1}], {A, i1}, {A, j1}] +  op[a, S[{B, i1}, {B, j1}], {B, i1}, {B, j1}]]
";

checkFields::usage="checkFields[expr] checks if all fields in expr are defined.\n

Example:\n
checkFields[op[S[{A, i1}, {B, j1}], {A, i1}, {B, j1}]]
";

checkIndices::usage="checkIndices[expr] checks if an index appears more often than twice in expr.\n

Example:
checkIndices[op[S[{A, i1}, {B, i1}], {A, i1}, {B, j1}]]
";

checkSyntax::usage="checkSyntax[expr] checks if expr has the correct syntax, i.e., op functions only contain propagators, vertices, fields, composite operators and regulator insertions and these quantities also have the correct arguments.\n

Examples:
checkSyntax[op[a,S[{A, i1}, {B, i2}], {A, i1}, {B, j1}]]

checkSyntax[dR[{A, i}, {A, j}, {A, l}]]
";

CO::usage="CO[{fieldCO, indexCO}, {field1, index1}, {field2, index2}, ...] represents a composite operator fieldCO with index indexCO of the fields fieldi with their indices indexi  in its symbolic form.
CO[fieldCO[mom, index1, index2, ...], fielda[momentuma, indexa1, indexa2, ...], fieldb[momentumb, indexb1, indexb2, ...], ..., explicit->True] represents a composite operator fieldCO with momentum mom and explicit indices indexi of the fields fieldji with their momenta momentumji and explicit indices indexjk in algebraic form.
The option explicit can have an arbitrary value.

Symbolic example: Composite operator G=phi^2
setFields[{G,phi}];
op[CO[{G,i}, {phi, i1}, {phi, i2}], {phi, i1}, {phi, i2}]/2

Algebraic example:
CO[oper[mom, inds], field1[mom1, inds1], field2[mom2, inds2], explicit -> True]

Example: Composite operator phi^a_i phi^a_j 
CO[G[p_, i_, j_], phi[p1_, is_, as_], phi[p2_, js_, bs_], explicit -> True]:=delta[as, bs] delta[i, is] delta[j, js]
";

complete::usage="complete is the default value for the option output of DSEPlot and RGEPlot.
";

complex::usage="complex is the field type of a complex bosonic field.
Properties of fields need to be set by setFields.
";

complexFieldQ::usage="complexFieldQ[expr] gives True if expr is a complex bosonic field.
";

connectedQ::usage="connectedQ[expr] gives True if expr is a connected diagram.";

COPlot::usage="COPlot[expr] plots expr which is an equation for correlation functions of composite operators.
COPlot[expr, styles] plots expr with formatting styles for the fields given in styles.

Possible options are:
	-) Options of Graph.
	
Example:
setFields[{phi, pp}];
action = {{phi, phi}};
Clear@F
F[j_] := Module[{i1, i2}, 
  op[CO[{pp, j}, {phi, i1}, {phi, i2}], {phi, i1}, {phi, i2}]/2]
FF = op[F[i], F[j]];
co = identifyGraphs[doCO[action, FF, onePIQ], {{pp, i}, {pp, j}}];
COPlot[co, {{phi, Black}}]
";

coSymbol::usage="coSymbol is an option of COPlot, DSEPlot, DSEPlotList and RGEPlot. It determines how to draw composite operators.

Possible values: boxSymbol, diskSymbol, triangleSymbol, diskTinySymbol, diskOpenSymbol, crossSymbol or a user-defined function which takes the coordinate of the regulator insertion as input.
Default value: triangleSymbol.

Example:
setFields[{A, FF}]; COPlot[
 op[CO[{FF, i}, {A, r1}, {A, s1}, {A, t1}, {A, u1}], 
  P[{A, r1}, {A, v1}], P[{A, s1}, {A, w1}], P[{A, t1}, {A, x1}], 
  P[{A, u1}, {A, y1}], 
  CO[{FF, j}, {A, v1}, {A, w1}, {A, x1}, {A, y1}]], {{A, Red}}, 
 coSymbol -> boxSymbol]";

countTerms::usage="countTerms[expr] counts the number of graphs appearing in expr.

Example:
countTerms[op[S[{A, i}, {A, j}]] +  1/2 op[S[{A, i}, {A, r1}, {A, s1}], V[{A, t1}, {A, u1}, {A, j}], P[{A, t1}, {A, r1}], P[{A, u1}, {A, s1}]]]
";

createDummyList::usage="createDummyList[n, dummyNames] creates ate least n names for dummy indices.
";

crossSymbol::usage="crossSymbol is a cross graphic used for bareVertexSymbol, coSymbol, regulatorSymbol or vertexSymbol.
";

deriv::usage="deriv[expr, {field, i}] differentiates expr with respect to {field, i}.
deriv[expr, fields] performs several derivatives with respect to the fields given in flis.
This function is used in doDSE.

Examples:
deriv[op[S[{A, r}, {A, s}, {A, t}], {A, r}, {A, s}, {A, t}], {A, i}]
deriv[op[S[{A, r}, {A, s}, {A, t}], {A, r}, {A, s}, {A, t}], {A, i},{A,j},{A,l}]
";

derivRGE::usage="derivRGE[expr, {field, i}] differentiates expr with respect to {field, i}.
deriv[expr, fields] perform several derivatives with respect to the fields given in flis.
This function is used in doRGE.

Examples:
derivRGE[op[V[{phi, i}, {phi, s}, {phi, t}], P[{phi, s}, {phi, t}]], {phi, j}]
derivRGE[op[V[{phi, i}, {phi, s}, {phi, t}], P[{phi, s}, {phi, t}]], {phi, j}, {phi, l}]
";

diagramTypes::usage="diagramTypes is an association containing all known diagram types.
A diagram type is defined by a list {n, {v1, v2, ...}}, where n is the loop number and the vi are the numbers of legs of all vertices.

Examples:
diagramTypes[{2, {4, 3, 3}}]
diagramTypes[{1, {3, 3, 3, 3}}]
";

disconnectedQ::usage="disconnectedQ[expr] gives True if expr is a disconnected diagram.
";

diskSymbol::usage="diskSymbol is a disk graphic used for bareVertexSymbol, coSymbol, regulatorSymbol or vertexSymbol.
";

diskOpenSymbol::usage="diskOpenSymbol is an open disk graphic used for bareVertexSymbol, coSymbol, regulatorSymbol or vertexSymbol.
";

diskTinySymbol::usage="diskTinySymbol is a tiny disk graphic used for bareVertexSymbol, coSymbol, regulatorSymbol or vertexSymbol.
";

doCO::usage = "doCO[ac, cf, [filter, opts]] derives the equation for the correlation function cf of composite operators with the action ac. filter are optional functions to select a subset of diagrams.

Example:
setFields[{phi, pp}];
action = {{phi, phi}};
Clear@F
F[j_] := Module[{i1, i2}, 
  op[CO[{pp, j}, {phi, i1}, {phi, i2}], {phi, i1}, {phi, i2}]/2]
FF = op[F[i], F[j]];
co = identifyGraphs[doCO[action, FF, onePIQ], {{pp, i}, {pp, j}}];
COPlot[co, {{phi, Black}}]
";

doDSE::usage="doDSE[ac, flis, [opts]] derives the DSE from the action ac for the fields contained in flis.
doDSE[ac, flis, props, [opts]] derives the DSE only with propagators contained in props given as {{field1a, field1b}, {field2a, field2b}.
doDSE[ac, flist, vtest, [opts]] derives the DSE only with vertices allowed by vtest.

Allowed propagators will be taken from ac if the props argument is not given. It is required, e.g., for actions with fields mixing at the two-point level but no corresponding tree-level propagator.

Examples:
Two-point DSE of an O(N) symmetric scalar theory
setFields[{phi}];
dse=doDSE[{{phi, phi}, {phi, phi, phi, phi}}, {phi, phi}]
DSEPlot[dse, {{phi,Black}}]

Two-point DSE of an O(N) symmetric scalar theory in the phase with broken symmetry. Two different ansaetze are given.
setFields[{phi}];
dse1 = doDSE[{{phi, phi}, {phi, phi, phi, phi}}, {phi, phi}, symmetry -> broken, ansatz :> {{phi, phi, phi, phi}}];
DSEPlot[dse1, {{phi, Black}}]
dse2 = doDSE[{{phi, phi}, {phi, phi, phi, phi}}, {phi, phi}, symmetry -> broken, ansatz :> {{phi, 6}}];
DSEPlot[dse2, {{phi, Black}}]

Three-gluon DSE of Landau gauge Yang-Mills theory with the dressed four-point vertices discarded by using a test function for all vertices
Clear@vTest; vTest[a_V] := Length@a < 4;
setFields[{A},{{c,cb}}];
dse = doDSE[{{A, A}, {c, cb}, {A, A, A}, {A, cb, c}, {A, A, A, A}}, {A, A, A}, {{A, A}, {c, cb}}, vTest];
DSEPlot[dse, {{A, Red}, {c, Green, Dashed}}]

Three-gluon DSE of Landau gauge Yang-Mills theory with the dressed four-point vertices discarded by specifying an ansatz for the effective action
setFields[{A},{{c,cb}}];
dse = doDSE[{{A, A}, {c, cb}, {A, A, A}, {A, cb, c}, {A, A, A, A}}, {A, A, A}, {{A, A}, {c, cb}}, ansatz -> {{A, A, A}, {A, cb, c}}];
DSEPlot[dse, {{A, Red}, {c, Green, Dashed}}]

Three-point DSE of a theory with bosonic fields A, phi, and phib which mix at the two-point level.
setFields[{A},{},{{phi,phib}}];
dse = doDSE[{{A, A}, {phi, phib}, {A, phi}, {A, phib}, {A, phib, phi}}, {A, A}, {{phi, phi}, {phib, phib}}]
DSEPlot[dse]
";

doGrassmannTest::usage="doGrassmannTest is an option of setSourcesZero. It ensures that the Grassmann number of each vertex is zero for each Grassmann field.
Only used for DSEs, as the ansatz for the action fully determines the allowed vertices for RGEs.
";

doRGE::usage="doRGE[ac, flis, [opts]] derives the RGE from the action ac for the fields contained in flis.
doRGE[ac, flis, props, [opts]] derives the RGE only with propagators contained in props given as {{field1a, field1b}, {field2a, field2b}.
doDSE[ac, flist, vtest, [opts]] derives the RGE only with vertices allowed by vtest.

Allowed propagators will be taken from ac if the props argument is not given. It is required, e.g., for actions with fields mixing at the two-point level and no corresponding tree-level propagator.

Examples:
Two-point RGE of an O(N) symmetric scalar theory in the symmetric phase
setFields[{phi}];
rge=doRGE[{{phi, phi}, {phi, phi, phi, phi}}, {phi, phi}]
RGEPlot[rge,{{phi,Black}}, output->forceEquation]

Two-point RGE of an O(N) symmetric scalar theory in the phase with broken symmetry
setFields[{phi}];
rge=doRGE[{{phi, phi}, {phi, phi, phi, phi}}, {phi, phi}, symmetry->broken]
RGEPlot[rge,{{phi,Black}}]

Three-gluon RGE of Landau gauge Yang-Mills theory with the four-point vertices discarded
Clear@vTest; vTest[a_V] := Length@a < 4;
setFields[{A},{{c,cb}}];
rge = doRGE[{{A, A}, {c, cb}, {A, A, A}, {A, cb, c}, {A, A, A, A}}, {A, A, A}, {{A, A}, {c, cb}}, vTest];
RGEPlot[rge, {{A, Red}, {c, Green, Dashed}}]

Three-point RGE of phi^3 theory with no regulator insertions
setFields[{phi}];
rge = doRGE[{{phi, phi}, {phi,  phi, phi}}, {phi, phi, phi}, tDerivative -> False]
RGEPlot[rge, {{phi, Black}}, output -> forceEquation]

Three-point RGE of a theory with bosonic fields A, phi, and phib which mix at the two-point level, i.e., additional propagators have to be given in an extra argument. No regulator insertions performed.
setFields[{A},{},{{phi,phib}}];
rge = doRGE[{{A, A}, {phi, phib}, {A, phi}, {A, phib}, {A, phib, phi}}, {A, A, A}, {{phi, phi}, {phib, phib}}, specificFieldDefinitions -> {A, phi, phib}, tDerivative -> False]
RGEPlot[rge]
";

dR::usage="dR[{field1, ind1}, {field2, ind2}] represents in symbolic form a regulator insertion, \[PartialD]_t R_k, where fieldi are fields and their indi indices.
dR[field1[mom1, inds1], field2[mom2, inds2], explicit -> True] represents a regulator insertion as needed by getAE. fieldi are fields, momi their momenta and indsi their full indices.

Example: Symbolic representation of a regulator insertion for gluons
dR[{A,i},{A,j}]

Example: Definition of regulator insertion for a scalar field with an O(N) index
dR[phi[p1,i], phi[p2,j], explicit -> True]:=delta[i,j] p1^2 dr[p1^2/k^2]
";

dummy::usage="dummy[i] represents a dummy index, i.e., an index over which is summed and integrated as appropriate. It is created by several functions of DoFun.
If the user requires dummy indices, the command insDummy[] should be used to guarantee the uniqueness of this variable.
";

DSEPlot::usage="DSEPlot[expr] plots a DSE expr.
DSEPlot[expr, fieldStyles] plots a DSE expr with the styles of the fields given by fieldStyles. The syntax is {{field1, style1}, {field2, style2}, ...}} where stylei are graphics primitives like colors suitable for Line.
DSEPlot[expr, n] or DSEPlot[expr, fieldStyles, n] plots the DSE expr with n graphs per row.
By default, blobs denote dressed quantities (with the exception of internal propagators), dots bare n-point functions and external fields are indicated by a circle.

Possible options are:
 -) Options of Graph.
 
Examples:
The gluon two-point DSE of Landau gauge Yang-Mills theory
setFields[{A},{{c,cb}}];
dse = doDSE[{{A, A}, {c, cb}, {A, cb, c}, {A, A, A}, {A, A, A, A}}, {A, A}];
DSEPlot[dse]

The gluon two-point DSE of Landau gauge Yang-Mills theory with gluons in red and ghosts dashed in green and four graphs per row
setFields[{A},{{c,cb}}];
dse = doDSE[{{A, A}, {c, cb}, {A, cb, c}, {A, A, A}, {A, A, A, A}}, {A, A}];
DSEPlot[dse,  {{A, Red}, {c, Green, Dashed}}, 4]

The graphs of the gluon two-point DSE of Landau gauge Yang-Mills theory in a list
setFields[{A},{{c,cb}}]; 
dse = doDSE[{{A, A}, {c, cb}, {A, cb, c}, {A, A, A}, {A, A, A, A}}, {A, A}];
DSEPlot[dse,  {{A, Red}, {c, Green}}, output -> List]

The complete two-point DSE of a free theory
setFields[{phi}];
dse = doDSE[{{phi, 2}}, {phi, phi}];
DSEPlot[dse,  {{phi, Black}}, output -> forceEquation]

Plotting a user-created expression with one external field
setFields[{phi}, {}, {}];
DSEPlot[op[S[{phi, i}, {phi, j}, {phi, l}, {phi, m}], P[{phi, l}, {phi, m}], {phi, j}], {{phi, Black}}]
";

DSEPlotList::usage="DSEPlotList[expr] returns a list of plots of each term in expr.
DSEPlot[expr, fieldStyles] returns a list of plots of each term in expr with the styles of the fields given by fieldStyles. The syntax is {{field1, style1}, {field2, style2}, ...}} where stylei are graphics primitives like colors suitable for Line.
DSEPlotList is the basic plotting function used by COPlot, DSEPlot, RGEPlot.
By default, blobs denote dressed quantities (with the exception of internal propagators), dots bare n-point functions, boxes regulator insertions, triangles composite operators and external fields are indicated by a circle.
By default, blobs denote dressed quantities (with the exception of internal propagators), dots bare n-point functions, boxes regulator insertions, triangles composite operators and external fields are indicated by a circle.
";

dummyCounter::usage="dummyCounter is a counter for dummy indices used by insDummy[].
";

dummyNames::usage="dummyNames is an option of createDummyList and sets the names for the dummy indices.
";

even::usage="even specifies that a field has only interactions with an even number of legs. See generateAction for details.
";

explicit::usage="explicit determines if the explicit form of a quantity is given. Propagators P, bare vertices S, dressed vertices V, composite operators CO and regulator insertions dR take it is an argument.
";

extractDiagramType::usage="extractDiagramType[diags, t] extracts diagrams of type t from diags. Known diagram types are stored in $diagramTypes.
";

factorStyle::usage="factorStyle is an option of COPlot, DSEPlot, DSEPlotList and RGEPlot. It determines the style of all text except indices and field labels.
Standard value: {FontSize:>16}.\n

Example:
setFields[{phi}];
rge = doRGE[{{phi, 4}}, {{phi, i}, { phi, j}}];
RGEPlot[rge, {{phi, Black}}, factorStyle -> {FontSize -> 20, Red, FontWeight -> Bold}]
";

fermion::usage="fermion is the field type of a Grassmann field.
Properties of fields need to be set by setFields.
";

fermionQ::usage="fermionQ[expr] returns True if expr is a Grassmann field.
";

fieldQ::usage="fieldQ[expr] yields True if expr is a field.
";

fieldType::usage="fieldType[expr] gives the field type of expr. Possible values are stored in $fieldTypes.
";

forceEquation::usage="forceEquation makes sure an equation is plotted for a single diagram.
";

generateAction::usage="generateAction[interacs] generates the action from interacs. Interactions are given as lists of the involved fields, e.g. {A,A,A}.
Symmetry factors are created automatically or can be given explicitly, e.g. {{A,A,A},6}.

The list of interactions can have the following elements:
 -) n-point functions as list of fields, e.g., {phi, phi} or {cb, c, A}
 -) A bosonic field and its maximal multiplicity, e.g., {phi, 4} will give two-, three- and four-point interactions.
 -) A bosonic field, its maximal multiplicity and the argument even to indicate that only interactions with an even number of fields involved should be taken into account, e.g., {phi, 4, even} will give two- and four-point interactions.
 -) A pair of bosonic complex fields or a pair of Grassmann fields and the maximal multiplicity of the pairs, e.g., {psi, psib, 2} will give the two- and the four-point functions.
 
Examples:
setFields[{A,phi},{{psi,psib}}];
generateAction[{{A,A},{A,A,A}}]
generateAction[{{phi, 4}}]
generateAction[{{phi, 4, even}}]
generateAction[{{psi, psib, 2}}]
generateAction[{{psi, psib}, {psib, psib, psi, psi}}]
";

get1PI::usage="get1PI[expr] returns the 1PI diagrams from expr.
";

getConnected::usage="getConnected[expr] returns the connected diagrams from expr.
";

getDiagramType::usage="getDiagramType[expr] returns the type of the diagram expr. Known diagram types are stored in $diagramTypes.
";

getDisconnected::usage="getDisconnected[expr] returns the disconnected diagrams from expr.
";

getInteractionList::usage="getInteractionList[ac] generates the list of interactions from a given symbolic action ac.\n

Example:
setFields[{A}];
getInteractionList[1/2 op[S[{A, r1}, {A, s1}], {A, r1}, {A, s1}] - 1/6 op[S[{A, u1}, {A, v1}, {A, w1}], {A, u1}, {A, v1}, {A, w1}]]
";
 
getLoopNumber::usage="getLoopNumber[expr] returns the number of loops of a diagram. If expr is a sum of diagrams, a list with the loop numbers is returned.

Example:
setFields[{phi}];
dse = doDSE[{{phi, 4}}, {phi, phi}];
getLoopNumber@dse
";

getNon1PI::usage="getNon1PI[expr] return the diagrams from expr which are not 1PI.
";

getSigns::usage="getSigns[expr] make signs from the auxiliary function sf explicit.
";

getVertexNumbers::usage="getVertexNumbers[expr] returns the number of vertices in expr.
";

grassmannQ::usage="grassmannQ[expr] return True if expr is a Grassmann field.\n
";

groupDiagrams::usage="groupDiagrams[expr] returns a list of diagrams in expr grouped by their types. Known diagram types are stored in $diagramTypes.";

identify::usage="identify is an option of doDSE and doRGE. It allows switching off the identification of diagrams.
";

identifyGraphs::usage="identifyGraphs[expr, extFields] adds up equivalent diagrams in expr. extFields are the external fields.
Note that this may not work for more than two loops.

Example:
setFields[{A}, {}, {}];
identifyGraphs[op[V[{A, i}, {A, r}, {A, s}, {A, j}], P[{A, r}, {A, s}]] + op[V[{A, i}, {A, j}, {A, s}, {A, t}], P[{A, s}, {A, t}]], {{A, i}, {A, j}}]
";

indexStyle::usage="indexStyle is an option for COPlot, DSEPlot, DSEPlotList and RGEPlot. It determines the style of all indices.
Standard value: {FontSize:>14}.

Example:
dse = doDSE[{{phi,phi}, {phi,phi,phi,phi}}, {phi, phi}];
DSEPlot[dse, {{phi, Black}}, indexStyle -> {FontSize -> 20, Blue, FontSlant -> Italic}]
";

insDummy::usage="insDummy[] returns unique dummy variable.

Example: Write down a graph using unique dummy variables
Module[{ind1=insDummy[],ind2=insDummy[]}, op[S[{phi,i},{phi,j},{phi,ind1},{phi,ind2}], P[{phi,ind1},{phi,ind2}]]]
";

intact::usage="intact is a value for the option symmetry.
"

odd::usage="odd specifies that a field can have interactions with an odd number of legs. Default value. See generateAction for details.
";

onePIQ::usage="onePIQ[expr] returns True if expr is a 1PI diagram.
";

op::usage="op[expr] can be used for symbolic and algebraic expressions representing combinations of propagators, vertices and so on.

Symbolic form:
Operator comprising (bare) vertices, propagators, composite operators, external fields and regulator insertions.
expr can be fields (denoted by {field, index}), bare vertices S, dressed vertices V, propagators P, regulator insertions dR or composite operators CO.
Summation and integration over mutliple indices is understood.
op automatically does some simplifications (see examples).
op is used in this way in the derivation of functional equations and when plotting them.

Examples:
op[S[{A,i},{A,r},{A,u}], P[{A,r},{A,s}], P[{A,u},{A,v}], V[{A,j}, {A,s}, {A,v}]]
op[0, V[{A, i}, {A, r}, {A, u}]]
op[2 S[{A, i}, {A, r}, {A, u}]]
op[V[{A, i}, {A, r}, {A, u}] - S[{A, i}, {A, r}, {A, u}], P[{A, r}, {A, u}]]

Algebraic form:
Operator comprising fields in the definition of physical actions.
The fields are given with all their indices in the form field[momentum, index1, index2, ...].
Summation and integration over mutliple indices is understood.
op is used in this way in getFR and convertAction.

Example: The two-point part of an O(N) symmetric scalar theory
setFields[{phi}];
convertAction[1/2 p^2 op[phi[p, i], phi[-p,i]]] 
";

orderFermions::usage="orderFermions[expr] orders Grassmann fields in expr canonically.

Note: orderFermions is deprecated and superseded by sortCanonical.

The canonical order is the following:
 -) vertices (V,S), regulator insertions (dR): antiFermions left of fermions
 -) propagators (P): antiFermions right (!) of fermions (In propagators the meaning of fermions and antiFermions is reversed for easier reading!)
 
Example:
setFields[{A}, {{c, cb}}, {}];
orderFermions[op[V[{c, i}, {cb, j}, {A, l}]]]
";

output::usage="output  is an option of COPlot, DSEPlot, DSEPlotList and RGEPlot. It determines in what form the output is given.

Possible values are:
 -) List: Gives a list of all graphs.
 -) forceEquation: Output in form of an equation, even if a single graph is plotted.
 -) complete (default): Output for several graphs in form of an equation and for a single graph as such.
 
Examples:
The graphs of the gluon two-point DSE of Landau gauge Yang-Mills theory in a list
setFields[{A},{{c,cb}}]; 
dse = doDSE[{{A, A}, {c, cb}, {A, cb, c}, {A, A, A}, {A, A, A, A}}, {A, A}];
DSEPlot[dse,  {{A, Red}, {c, Green}}, output -> List]

The complete two-point RGE of an O(N) symmetric scalar theory in the symmetric phase
setFields[{phi}];
dse = doRGE[{{phi, phi}, {phi,phi,phi,phi}}, {phi, phi}];
RGEPlot[dse,  {{phi, Black}}, output -> forceEquation]
";

P::usage="P[{field1, index1}, {field2, index2}] represents a dressed propagator of the fields fieldi with their indices indexi in its symbolic form.
P[field1[momentum1, index1a, index1b, ...], field2[momentum2, index2a, index2b, ...], explicit->True] represents a dressed propagator of the fields fieldi with their momenta momentumi and explicit indices indexij in algebraic form.
The option explicit can have an arbitrary value.

Symbolic example:
P[{A,i},{A,j}]

Algebraic example:
P[field1[mom1, inds1], field2[mom2, inds2], explicit -> True]

Example: Definition of a dressed propagator for a scalar field with an O(N) index
P[phi[p1_,i_], phi[p2_,j_], explicit -> True]:=delta[i,j] D[p1^2]/p1^2
";

propagatorCreationRules::usage="propagatorCreationRules is an option of setSourcesZero. It is used to distinguish between DSEs and RGEs.
";

regulatorBox::usage="regulatorBox is superseded by boxSymbol.
";

regulatorCross::usage="regulatorCross is superseded by crossSymbol.
";

regulatorSymbol::usage="regulatorSymbol  is an option of COPlot, DSEPlot, DSEPlotList and RGEPlot. It determines how to draw regulator insertions.

Possible values: boxSymbol, diskSymbol, triangleSymbol, diskTinySymbol, diskOpenSymbol, crossSymbol or a user-defined function which takes the coordinate of the regulator insertion as input.
Default value: boxSymbol.

Example:
setFields[{phi}, {}, {}];
RGEPlot[1/2 op[dR[{phi, r1}, {phi, s1}], P[{phi, t1}, {phi, r1}], P[{phi, s1}, {phi, v1}], V[{phi, i}, {phi, j}, {phi, v1}, {phi, t1}]], {{phi, Black}}, regulatorSymbol -> ({Text[\"Here comes the regulator.\", #]} &)]
";

replacedField::usage="replacedField[expr] represents a field expr replaced by the field and the propagator with a derivative. Only an intermediate dummy.
";

replaceFields::usage="replaceFields[expr] replaces fields in expr by the field and the propagator with a derivative.

Example:
setFields[{A}];
replaceFields[op[S[{A, i}, {A, r}, {A, s}], {A, r}, {A, s}]]
";

resetDummy::usage="resetDummy[] resets the counter in the dummy function.
Note: Should be used with care, because the uniqueness of dummy indices is not guaranteed after using resetDummy.

Example:
{insDummy[], insDummy[], resetDummy[], insDummy[]}
";

RGEPlot::usage="RGEPlot[expr] plots an RGE expr.
RGEPlot[expr, fieldStyle] plots an RGE expr with  the styles of the fields given by fieldStyles. The syntax is {{field1, style1}, {field2, style2}, ...}} where stylei are graphics primitives like colors suitable for Line.
RGEPlot[expr, n] or RGEPlot[expr, fieldStyles, n] plots the RGE expr with n graphs per row.
By default, blobs denote dressed quantities (with the exception of internal propagators), dots bare n-point functions, boxes regulator insertions and external fields are indicated by a circle.

Possible options are:
 -) Options of Graph.
 
Examples:
The gluon two-point RGE of Landau gauge Yang-Mills theory
setFields[{A},{{c,cb}}];
dse = doRGE[{{A, A}, {c, cb}, {A, cb, c}, {A, A, A}, {A, A, A, A}}, {A, A}];
RGEPlot[dse]

The gluon three-point RGE of Landau gauge Yang-Mills theory with gluons in red and ghosts dashed in green and four graphs per row
setFields[{A},{{c,cb}}];
dse = doRGE[{{A, A}, {c, cb}, {A, cb, c}, {A, A, A}, {A, A, A, A}}, {A, A}];
RGEPlot[dse,  {{A, Red}, {c, Green, Dashed}}, 4]

The graphs of the gluon two-point RGE of Landau gauge Yang-Mills theory in a list
setFields[{A},{{c,cb}}]; 
dse = doRGE[{{A, A}, {c, cb}, {A, cb, c}, {A, A, A}, {A, A, A, A}}, {A, A}];
RGEPlot[dse,  {{A, Red}, {c, Green}}, output -> List]

The complete two-point RGE of an O(N) symmetric scalar theory in the symmetric phase
setFields[{phi}];
dse = doRGE[{{phi, 100, even}}, {phi, phi}];
RGEPlot[dse,  {{phi, Black}}, output -> forceEquation]

Plotting a user-created expression with one external field
setFields[{phi}, {}, {}];
RGEPlot[op[S[{phi, i}, {phi, j}, {phi, l}, {phi, m}], P[{phi, l}, {phi, m}], {phi, j}], {{phi, Black}}]
";

S::usage="S[{field1, index1}, {field2, index2}, ...] represents a bare vertex, i.e., an expansion coefficient of the action, of the fields fieldi with their indices indexi in its symbolic form.
S[field1[momentum1, index1a, index1b, ...], field2[momentum2, index2a, index2b, ...], ..., explicit->True] represents a bare vertex of the fields fieldi with their momenta momentumi and explicit indices indexij in algebraic form.
The option explicit can have an arbitrary value.

Symbolic example:
S[{A,i},{A,j},{A,k}]

Algebraic example:
S[field1[mom1, inds1], field2[mom2, inds2], explicit -> True]

Example: Definition of a bare vertex for a scalar field with an O(N) index
S[phi[p1_,i_], phi[p2_,j_], phi[p3_,l_], phi[p4_,m_], explicit -> True]:=g (delta[i,j]delta[l,m]+delta[i,l]delta[j,m]+delta[i,m]delta[j,l])
";

sE::usage:="sE is identical to shortExpression."

setFields::usage = "setFields[bos, ferm, comp] sets the properties of the real bosonic fields bos, the Grassmann fields ferm and the complex fields comp.
setFields[bos, ferm] sets the properties of the real bosonic fields bos and the Grassmann fields ferm.
setFields[bos] sets the properties of the real bosonic fields bos.
The real bosonic fields are given as lists, the Grassmann and complex fields as pairs of lists. 

Example: Definition of a bosonic field A, a pair of anti-commuting fields c and cb and a pair of bosonic complex fields phi and phib.
setFields[{A}, {{c, cb}}, {{phi, phib}}];
bosonQ /@ {A, c, cb, phi, phib}
fermionQ /@ {A, c, cb, phi, phib}
antiComplexFieldQ /@ {A, c, cb, phi, phib}
antiField/@{A,c,cb,phi,phib}
";

setSourcesZero::usage="setSourcesZero[expr, ac, extLegs] sets the sources in expr with external legs extLegs to zero, i.e., only physical propagators and vertices for the action ac are left.
setSourcesZero[expr, ac, extLegs, ownAllowedPropagators] sets the sources in expr with external legs extLegs to zero with ownAllowedPropagators a list of propagators allowed. Given in the form {{field1a, field1b}, {field2a, field2b}, ...}.
setSourcesZero[expr, ac, legs, ownAllowedPropagators, vertexTest, opts] sets the sources in expr with external legs extLegs to zero with vertexTest a function to determine if a vertex should be kept.
This function is not for RGEs.

Examples:
One external field
setSourcesZero[op[S[{A, i}, {A, j}, {A, r}], {A, r}], {{A, A}, {A, A, A}}, {{A, A}}]

Replace dummy fields by physical fields
setSourcesZero[op[S[{A, i}, {A, r}, {A, s}, {A, t}], P[{A, r}, {$dummyField, u}], P[{A, s}, {$dummyField, v}], P[{A, t}, {$dummyField, w}], V[{$dummyField, u}, {$dummyField, v}, {$dummyField, w}, {A, j}]], {{A, A}, {A, A, A}},{{A, A}}]

Replace dummy fields by real fields and apply a test for the resulting vertices
Clear@vTest; vTest[a_V] := Length@a < 4;
setSourcesZero[op[S[{A, i}, {A, r}, {A, s}, {A, t}], P[{A, r}, {$dummyField, u}], P[{A, s}, {$dummyField, v}], P[{A, t}, {$dummyField, w}], V[{$dummyField, u}, {$dummyField, v}, {$dummyField, w}, {A, j}]], {{A, A}, {A, A, A}},{{A, A}}, vTest]
";
  
setSourcesZeroRGE::usage="setSourcesZeroRGE[expr, ac, extLegs] sets the sources in expr with external legs extLegs to zero, i.e., only physical propagators and vertices for the action ac are left.
setsetSourcesZeroRGEpr, ac, extLegs, ownAllowedPropagators] sets the sources in expr with external legs extLegs to zero with ownAllowedPropagators a list of propagators allowed. Given in the form {{field1a, field1b}, {field2a, field2b}, ...}.
setsetSourcesZeroRGEpr, ac, legs, ownAllowedPropagators, vertexTest, opts] sets the sources in expr with external legs extLegs to zero with vertexTest a function to determine if a vertex should be kept.
This function is for RGEs.

Examples:
One external field
setSourcesZeroRGE[op[V[{A, i}, {A, j}, {A, r}], {A, r}], {{A, A}, {A, A, A, A}}, {{A, A}}]

Replace dummy fields by physical fields
setSourcesZeroRGE[op[V[{A, i}, {A, r}, {A, s}, {A, t}], P[{A, r}, {$dummyField, u}], P[{A, s}, {$dummyField, v}], P[{A, t}, {$dummyField, traceIndex2}], V[{$dummyField, u}, {$dummyField, v}, {$dummyField, traceIndex2}, {A, j}]], {{A, A}, {A, A, A, A}},{{A, A}}]

Replace dummy fields by real fields and apply a test for the resulting vertices
Clear@vTest; vTest[a_V] := Length@a < 4;
setSourcesZeroRGE[op[V[{A, i}, {A, r}, {A, s}, {A, t}], P[{A, r}, {$dummyField, u}], P[{A, s}, {$dummyField, v}], P[{A, t}, {$dummyField, traceIndex2}], V[{$dummyField, u}, {$dummyField, v}, {$dummyField, traceIndex2}, {A, j}]], {{A, A}, {A, A, A}},{{A, A}}, vTest]
";

sf::usage="sf[field1, {field2, field3, ...}] encodes the sign for Grassmann fields. It is -1 if field1 is a Grassmann field and there is an odd number of Grassmann fields in the second argument.
Fields are given as {field, index}.
Some simplifications are done automatically. 
Signs are made explicit with getSigns.

Examples:
setFields[{A}, {{c, cb}}];
getSigns[sf[{cb, i}, {{cb, j}}]]
getSigns[sf[{cb, i}, {{cb, j}, {c, k}}]]
sf[{A, i}, {{cb, j}}]
";

shortExpression::usage="shortExpression[expr] writes a symbolic expression expr containing propagators, vertices and so on into a shorter form using $bareVertexSymbol, $vertexSymbol, $regulatorInsertionSymbol, $compOpSymbol and $propagatorSymbol for representation.
shortExpression[exp, opts] writes expr with the formatting options given in opts.
The function sE is identical to shortExpression.

Example:
shortExpression[1/2 op[S[{A, i}, {A, r}, {A, s}], V[{A, t}, {A, u}, {A, j}], P[{A, t}, {A, r}], P[{A, u}, {A, s}]], Red, FontSize -> 20]\n
";

sortCanonical::usage="sortCanonical[expr] orders the fields in vertices and propagators in expr in a canonical way:
-) Anti-fermions left of fermions in vertices.
-) Fermions left of anti-fermions in propagators. This is due to the definition of propagators, which show the anti-fields instead of the fields to allow easier identification with the corresponding vertex legs.
-) External fields fields ordered by list of derivatives.
-) Internal fields ordered by connection to external fields.

Example:
setFields[{A}, {{c, cb}}];
sortCanonical[op[V[{c, a}, {cb, b}, {A, i}], P[{cb, b}, {c, f}], 
  P[{c, a}, {cb, g}], V[{cb, g}, {c, f}, {A, j}]], {{A, i}, {A, j}}]
";

S::usage="S[{field1, index1}, {field2, index2}, ...] represents a bare vertex, i.e., an expansion coefficient of the action, of the fields fieldi with their indices indexi in its symbolic form.
S[field1[momentum1, index1a, index1b, ...], field2[momentum2, index2a, index2b, ...], ..., explicit->True] represents a bare vertex of the fields fieldi with their momenta momentumi and explicit indices indexij in algebraic form.
The option explicit can have an arbitrary value.

Symbolic example:
S[{A,i},{A,j},{A,k}]

Algebraic example:
S[field1[mom1, inds1], field2[mom2, inds2], explicit -> True]

Example: Definition of a bare vertex for a scalar field with an O(N) index
S[phi[p1_,i_], phi[p2_,j_], phi[p3_,l_], phi[p4_,m_], explicit -> True]:=g (delta[i,j]delta[l,m]+delta[i,l]delta[j,m]+delta[i,m]delta[j,l])
";

sortDummies::usage="sortDummies[expr] replaces the dummy indices by shorter dummies making the expression thus easier to read.
This function is automatically applied by some functions.

Example:
sortDummies[op[S[{phi, i100}, {phi, j1}, {phi, myInternalIndexWithALongName}, {phi, myExternalIndexWithALongNames}], P[{phi, i100}, {phi, j1}], {phi, myInternalIndexWithALongName}]]
";

sourcesZero::usage="sourcesZero is an option of doCO, doDSE and doRGE. It allows avoiding setting the sources to zero. 
";

specificFieldDefinitions::usage="specificFieldDefinitions was removed in DoFun 3.0.0.
";

symmetry::usage="symmetry is an option of doCO, doDSE and doRGE. It determines if there is broken symmetry in the theory.
Option of doDSE and doRGE.

Possible values:
 -) broken
 -) intact (default)
";

superField::usage="superField is the field type of a super field representing all possible fields.
";

tDerivative::usage="tDerivative is an option of doRGE. It determines if the regulator insertion is performed.
Default: True.
";

traceIndex1::usage="traceIndex1 is a dummy index used by doRGE.
";

traceIndex2::usage="traceIndex2 is a dummy index used by doRGE.
";

traceIndices::usage="traceIndices are dummy indices used by doRGE.
";

triangleSymbol::usage="triangleSymbol is a triangle graphic used for bareVertexSymbol, coSymbol, regulatorSymbol or vertexSymbol.
";

type::usage="types is an option of COPlot, DSEPlot, DSEPlotList and RGEPlot. It serves to tell the underlying plot function DSEPlotList how to plot the left-hand side of the equation.

Possible values are:
	-) \"CO\"
	-) \"DSE\"
	-) \"RGE\"
";

userEvenFields::usage="userEvenFields is an option of doDSE and doRGE. It allows providing fields which can only appear in even numbers in vertices.
";

V::usage="V[{field1, index1}, {field2, index2}, ...] represents a dressed vertex of the fields fieldi with their indices indexi  in its symbolic form.
V[field1[momentum1, index1a, index1b, ...], field2[momentum2, index2a, index2b, ...], ..., explicit->True] represents a bare vertex of the fields fieldi with their momenta momentumi and explicit indices indexij in algebraic form.
The option explicit can have an arbitrary value.

Symbolic example:
V[{A,i},{A,j},{A,k}]

Algebraic example:
V[field1[mom1, inds1], field2[mom2, inds2], explicit -> True]

Example: Definition of a dressed vertex for a scalar field with an O(N) index
V[phi[p1_,i_], phi[p2_,j_], phi[p3_,l-], phi[p4_,m_], explicit -> True]:=g (delta[i,j]delta[l,m]+delta[i,l]delta[j,m]+delta[i,m]delta[j,l])
";

vertexSymbol::usage="vertexSymbol is an option of COPlot, DSEPlot, DSEPlotList and RGEPlot. It determines how to draw dressed vertices.

Possible values: boxSymbol, diskSymbol, triangleSymbol, diskTinySymbol, diskOpenSymbol, crossSymbol or a user-defined function which takes the coordinate of the regulator insertion as input.
Default value: diskSymbol.

Example:
setFields[{phi}, {}, {}];
op[P[{phi, t1}, {phi, v1}], V[{phi, i}, {phi, j}, {phi, v1}, {phi, t1}]], {{phi, Black}}, vertexSymbol -> triangleSymbol]
";





(* ::Section:: *)
(* Options *)


(* define the names for the dummy indices; later numbers will be added *)
Options[createDummyList]={dummyNames->{Global`r,Global`s,Global`t,Global`u,Global`v,Global`w,Global`x,Global`y,Global`z}};

Options[DSEPlot]={
	bareVertexSymbol->diskTinySymbol,
	coSymbol->triangleSymbol,
	factorStyle->{FontSize:>16},
	indexStyle->{FontSize:>14},
	output->complete,
	type->"DSE",
	regulatorSymbol->boxSymbol,
	vertexSymbol->diskSymbol,
	(* options of Graph *)
	ImageSize->100};

Options[DSEPlotList]:=Options[DSEPlot];

(* since RGEPlot and COPlot rely on DSEPlotList, the options should be the same; changing the former options may be ignored by some function, since they take the options of DSEPlot *)
Options[RGEPlot]:=Join[{type->"RGE"}, DeleteCases[Options[DSEPlot], type->_]];

Options[COPlot]:=Join[{type->"CO"}, DeleteCases[Options[DSEPlot], type->_]];

Options[setSourcesZero]={doGrassmannTest->True, propagatorCreationRules->DSERules};

Options[doDSE]={sourcesZero:> True, identify:> True, ansatz->{}, userEvenFields->{}};

Options[doRGE]={sourcesZero:> True, identify:> True, tDerivative -> True, 
	symmetry->intact, userEvenFields->{}};

Options[shortExpression]={FontSize->16};




(* ::Section:: *)
(* Variables *)


(* the standard symbol for a bare vertex used in shortExpression *)
$bareVertexSymbol=S;

(* the standard symbol for a propagator used in shortExpression *)
$propagatorSymbol=\[CapitalDelta];

(* the standard symbol for a vertex used in shortExpression *)
$vertexSymbol=\[CapitalGamma];

(* the standard symbol for a composite operator used in shortExpression *)
$compOpSymbol=\[Omicron];

(* the standard symbol for a regulator insertion used in shortExpression *)
$regulatorInsertionSymbol=R;
(* list of indices automatically used for external vertices *)


$externalIndices={Global`i1,Global`i2,Global`i3,Global`i4,Global`i5,Global`i6,Global`i7,Global`i8,Global`i9,Global`i10};
Protect/@$externalIndices;

(* list of all field types *)
$fieldTypes={boson, fermion, antiFermion, complex, antiComplex, superField};

(* the standard super field *)
$dummyField = \[Phi]; 
Evaluate@$dummyField /: fieldType[$dummyField] := superField;


(* fermionic dummy fields *)
(* removed after version 3.0.1 *)




(* Define various diagram types. 
The first entry specifies the number of loops, the second gives an ordered list of vertex multiplicities.
For example, {1, {3, 4}} represents a one-loop diagram with two three-point functions. *)
diagramTypes = <|(* two-point *)
   {1, {3, 3}} -> "oneLoop",
   {1, {4}} -> "tadpole",
   {2, {4, 3, 3}} -> "sunset",
   {2, {4, 4}} -> "squint",
   (* three-point *)
   {1, {3, 3, 3}} -> "triangle",
   {1, {3, 3, 3}} -> "triangle3",
   {1, {3, 4}} -> "swordfish",
   {1, {3, 4}} -> "swordfish3",
   (* four-point *)
   {1, {3, 3, 3, 3}} -> "box",
   {1, {3, 3, 4}} -> "triangle4",
   {1, {4, 4}} -> "swordfish4",
   {1, {3, 5}} -> "fivePoint4"|>;
   
(* Known classes of diagrams *)
$diagramTypes = Values@diagramTypes;

(* abbreviationa for some functions *)
sE=shortExpression;


(* switch for sign convention;
when $signConvention=1, the DSE are derived with the convention that 1PI vertex functions (n>2)
are the negative derivative of the effective action;
if it is -1, then the "mathematical" definition applies that the 1PI vertex functions
are the positive derivatives *)
$signConvention=-1;


(* open indices that are closed by the trace *)
traceIndices={traceIndex1,traceIndex2};




(* ::Section:: *)
(* Start private *)


Begin["`Private`"];

(*
private functions (alphabetic)
	arrowLine
	createDummyList
	createPropagatorRules
	derivAll
	derivAllRGE
	derivField
	derivPropagator
	derivPropagators
	derivPropagatorsdt
	derivVertex
	DSEPlotCompare
	DSEPlotGrid
	DSEPlotList
	evenBosonTest
	firstDerivReplacement
	getDirectedFieldsList
	getGraphCharacteristic
	getNeighbours
	getSignature
	grassmannTest
	indicesTest
	indicesTwice
	innerFermionsCanonicalQ
	insertRegulator
	newVertex
	newVertexSum
	opTest
	plugInFieldsV
	propagatorTest
	regulatorInsertionTest
	replacementCalc
	replacementCalcStep
	replaceVertices
	replaceVerticesBrokenPhase
	RGEPlotGrid
	setSourcesZeroDo
	shortExpressionSingle
	shortExpressionStyle
	vertexDummies
	vertexTest

*)




(* ::Section:: *)
(* Messages *)


checkAll::syntax="There was a syntax error in checkAll.\n
Make sure the input has the form of ... .\n
The expression causing the error is `1`.";


checkFields::ok="All fields are defined correctly.";

checkFields::undefinedField="The expression(s) in `1` is/are not defined as field(s). Use setFields or generateAction to do so.";


checkIndices::ok="No indices appear more often than twice.";

checkIndices::multipleIndices="The index `2` appears more than twice in `1`.";



checkAction::error="There is an error in `1`: The indices `2` do not appear twice. In the action all indices should be summed over.";

checkAction::ok="All summations ok.";



checkSyntax::op="There is a syntax error in `1`.";

checkSyntax::propagator="There is a syntax error in the propagator `1`.";

checkSyntax::vertex="There is a syntax error in the vertex `1`.";

checkSyntax::regulatorInsertion="There is a syntax error in the regulator insertion `1`.";

checkSyntax::ok="The syntax seems to be ok.";


countTerms::syntax="There was a syntax error in countTerms.\n
Make sure the input has the form of countTerms[expr_]. For more details use ?countTerms.\n
The expression causing the error is `1`.";


setFields::dummyField="Fields must not contain the dummy field "<>ToString@$dummyField<>". Rename your field.";

setFields::syntax="There was a syntax error in setFields.\n
Make sure the input has the form of setFields[bosons_List, [fermions_List, complexFields_List]]. For more details use ?setFields.\n
The expression causing the error is `1`.";


deriv::syntax="There was a syntax error in deriv.\n
Make sure the input has the form of deriv[expr_, flis_]. For more details use ?deriv.\n
The expression causing the error is `1`.";


derivAll::index="Use another index for differentiation. `1` is already contained in the action.";


DSEPlotList::syntax="There was a syntax error in DSEPlotList.\n
Make sure the input has the form of DSEPlotList[expr_,flis_List[,plotRules_],opts___]. For more details use ?DSEPlotList.\n
The expression causing the error is `1`.";

DSEPlot::fieldsUndefined="There appear to be undefined fields in the expression you want to plot.\n
Define the fields with the function setFields or rederive the expression with doRGE or doDSE, respectively.
The definition is then done automatically.\n
The expression causing the error is `1`";

DSEPlot::syntax="There was a syntax error in DSEPlot.\n
Make sure the input has the form of DSEPlot[expr_,flis_List[,plotRules_],opts___]. For more details use ?DSEPlot.\n
The expression causing the error is `1`.";

DSEPlot::plotRules="The style definitions for the fields are not properly defined.
They have to have the form {{field1, style1}, {field2, style2}, ...}.
The faulty style definitions are `1`";

RGEPlot::syntax="There was a syntax error in RGEPlot.\n
Make sure the input has the form of RGEPlot[expr_,flis_List[,plotRules_],opts___]. For more details use ?RGEPlot.\n
The expression causing the error is `1`.";

RGEPlot::plotRules=DSEPlot::plotRules;

generateAction::syntax="There was a syntax error in generateAction.\n
Make sure the input has the form of generateAction[interactions_List[,fields_List]]. For more details use ?generateAction.\n
The expression causing the error is `1`.";

generateAction::fieldsNotSet="There was a syntax error in generateAction.\n
The argument of generateAction contains expressions which are not fields. Use setFields to define them as fields. For more details use ?generateAction.\n
The expression causing the error is `1`.";

getLoopNumber::syntax="There was a syntax error in getLoopNumber.\n
Make sure the input has the form of getLoopNumber[expr_].\n
The expression causing the error is `1`.";

getInteractionList::syntax="There was a syntax error in getInteractionList.\n
Make sure the input has the form of getInteractionList[action_]. For more details use ?getInteractionList.\n
The expression causing the error is `1`.";

identifyGraphs::syntax="There was a syntax error in identifyGraphs.\n
Make sure the input has the correct form. For more details use ?identifyGraphs.\n
The expression causing the error is `1`.";

setFields::noList="There is a syntax error in setFields, because where a pair of Grassmann or complex bosonic fields was expected, something else was encountered.\n
Make sure to group fermions and bosonic complex fields in pairs, e.g., {{psi, psib}, {c,cb}}."

doDSE::syntax="There was a syntax error in doDSE.\n
Make sure the input has the form of doDSE[L_, flis_List, [propagators_List,] vertexFunction_Symbol]. For more details use ?doDSE.\n
The expression causing the error is `1`.";

doDSE::fieldAsIndex="Do not use field names as indices. Change the index/indices `1`.";

doDSE::tooFewExternalIndices="The list $externalIndices=`1` does not have enough indices. It must have at the last as many as there are legs of the diagram. Either change external indices or add indices to the fields manually. See ?doDSE for details on how to do this.";

doDSE::noAnsatz="This is no error message, only a warning.
The option symmetry->broken was set, but no ansatz was given. Are you sure this was done on purpose?
Use the option ansatz in doDSE to specify which interactions should be allowed.";

doRGE::syntax="There was a syntax error in doRGE.\n
Make sure the input has the form of doRGE[...]. For more details use ?doRGE.\n
The expression causing the error is `1`.";

doRGE::noTruncation="For a derivation of an RGE in the symmetry broken phase a truncated average effective action is required as the calculation
cannot be done up to arbitrary order. Do this by restricting the action to the vertices contained in the truncation.";


orderFermions::superseded="Note: orderFermions is deprecated and superseded by sortCanonical. It is no longer updated.";


RGEPlot::fieldsUndefined=DSEPlot::fieldsUndefined;

setFields::noPair = 
  "The fermions and/or complex fields are not given as pairs. The \
input was `1` and `2`. Each should be of the form {{field1, \
anti-field1},...}.";

setFields::syntax = 
  "There was a problem with the syntax of setFields.";
  
setSourcesZero::syntax="There was a syntax error in setSourcesZero.\n
Make sure the input has the form of setSourcesZero[expr_ [,propagators_List], vertexTest_Symbol]. Fore more details use ?setSourcesZero.\n
The expression causing the error is `1`.";

shortExpression::syntax="There was a syntax error in shortExpression.\n
Make sure the input has the form of shortExpression[expr_, opts___]. For more details use ?shortExpression.\n
The expression causing the error is `1`.";

sE::syntax=shortExpression::syntax;

sortCanonical::dummyFields="Warning: Expression was not sorted by sortCanonical, because it contains dummy fields. Make sure the sources are set to zero.
The expression is `1`";

DSEPlotList::multiPropagators="Note: More than two different propagators connecting the same vertices cannot be displayed properly. The plot might contain wrong styles for these propagators.\n";




(* ::Section:: *)
(* General Functions *)

(* op cleared for DoFun3 *)
(* the op functions hold together vertices and propagators; it stands for summations/integration over indices *)
op[a___,b_ ?NumericQ c_,d___]:=b op[a,c,d];

op[a___,Plus[b_, c_],d___]:=op[a,b,d]+op[a,c,d];

op[a___,b_?NumericQ,c___]:=b op[a,c];

(* property similar to Flat, but without side effects on pattern matching *)
op[a___,op[b___], c___]:=op[a,b,c];


(* dummies cleared for DoFun3 *)
(* Dummies: They are used to represent contracted indices.
They are realized by the function dummy[x], where x is a running number.
When the package is loaded, x is set to 0.
It can be rest using resetDummy.
Dummies are inserted using insDummy[].
For nicer output, sortDummies[] can replace dummies by other dummy variables.
The standard choice is r1, r2, ..., z1, z2, ..., depending on how many are needed.
This list is created with createDummyList[]. *)

dummyCounter=0;


resetDummy[]:=(dummyCounter=0;);


insDummy[]:=(dummyCounter++;dummy[dummyCounter]);


createDummyList[a_Integer,type_,opts___]:=Module[{dNames},
	dNames=Join[{opts},type/.Options@createDummyList];
	ToExpression/@Flatten@Outer[
	StringJoin,ToString/@dNames,ToString/@Range[Ceiling[a/Length@dNames]]]
];


sortDummies[a_Plus]:=sortDummies/@a;

sortDummies[a_Times]:=sortDummies/@a;

sortDummies[a_List]:=sortDummies/@a;

sortDummies[a_?NumericQ]:=a;

sortDummies[a_op]:=Module[
{b,inds,exprDummies,newDummies,newDummyRules,dummyList},

(* discard sf terms, since their indices are not to be contracted *)
b = a/.sf[_,_]:>1;

(* get all indices from the epxression, do not consider indices in sf *)
(*inds=Transpose[b[[Sequence@@#]]&/@Position[b,{_,ind_}]][[2]];*)
inds = Cases[b/.sf[_,_]:>1, {c_?fieldQ, ind_} :> ind, Infinity];

(* determine the dummy indices; only take the first appearance, don't change order, i.e. don't use Union *)
(*exprDummies=
Cases[{#,Count[inds,#]}&/@inds,{g_,2}:> g]/.{f___,g_,c___,g_,e___}:> {f,g,c,e};*)
exprDummies = Cases[Tally[inds],{c_,2}:>c];

(* create the dummy list *)
dummyList=createDummyList[Length@exprDummies,dummyNames];

(* determine rules for replacement *)
newDummies=Take[dummyList,Length@exprDummies];
newDummyRules=Thread[Rule[exprDummies,newDummies]];

a/.newDummyRules
];



(* create from a list of interactions from the action in the usual form;
an interaction is given as a list of the involved fields *)

(* in case a "true" action is given, strip it down to its parts relevant for doDSE/doRGE *)
generateAction[action_,rest___]:=
	generateAction[(Cases[Replace[#, op[b__] :> {op[b]}(*in case a=op[___]*)], op[___], \[Infinity]] /. op :> Sequence)[[All, 0]] & /@ List@@action,rest];

(* properties of fields are not set yet *)
(*generateAction[interactions_List]:=Message[generateAction::fieldsNotSet,interactions];*)

generateAction[interactions_List]:=Module[
{interactions2, fields, autoList, userList, allList, factor},

(* in case the fields are given by their order and symmetry *)
interactions2=interactions/. {
	{Q_, order_Integer, even} :>Sequence@@NestList[Drop[#, 2] &, Table[Q, {order}], Floor[order/2] - 1],
	{Q_, order_Integer, odd}|{Q_, order_Integer} :>Sequence@@NestList[Drop[#, 1] &, Table[Q, {order}], order - 2 (* go only down to propagator *)],
	{Q_,Qb_, order_Integer} :>Sequence@@NestList[Drop[Drop[#, 1], -1] &, Flatten@Transpose[Reverse /@ Table[{Q, Qb}, {order}]], Floor[order] - 1](*,
	{Q_,even}:>Sequence@@NestList[Drop[#, 2] &, Table[Q, {100}], Floor[100/2] - 1]*)
	};

(* add the numerical factors and the dummy indices;
the user may give his own prefactors, so split in two parts: user defined and no user factors *)

userList=Cases[interactions2,{fields_List, factor_}];
autoList=Complement[interactions2, userList];

userList=userList/.field_?fieldQ :> {field, insDummy[]};

autoList = Function[ia, 
   Transpose@(Tally@ia) /. {localFields_List, 
      localFactors_List} :> {ia /. field_?fieldQ :> {field, insDummy[]}, 
      1/Times @@ Factorial /@ localFactors}] /@ autoList;

(*(Plus @@ (#[[2]] op[S[Sequence @@ #[[1]]], Sequence @@ Reverse@#[[1]]] & /@ Join[autoList, userList] // sortDummies))
  /.a_S?(Length@#>2&):>-$signConvention a*)
  (* join lists *)
  allList = Join[autoList, userList];
  
  (* generate action; for vertices change the order of the fields to account for the correct order of fermions;
  for propagators not necessary because they are defined with field/anti-fields exchanged *)
  (*(Plus @@ (#[[2]] op[S[Sequence @@ #[[1]]], Sequence @@ #[[1]]] & /@ allList // sortDummies))
  /.a_S?(Length@#>2&):>-$signConvention a /. op[b_S,c___List]:>op[b,Sequence@@Reverse[{c}]]/;Length[{c}]>2*)
(*  Print[(Plus @@ (#[[2]] op[S[Sequence @@ #[[1]]], Sequence @@ #[[1]]] & /@ allList // sortDummies))
  /.a_S?(Length@#f>2&):>-$signConvention a ];*)
  (Plus @@ (#[[2]] op[S[Sequence @@ #[[1]]], Sequence @@ #[[1]]] & /@ allList // sortDummies))
  /.a_S?(Length@#>2&):>-$signConvention a /. op[b_S,c___List]:>op[b,Sequence@@Reverse[{c}]]/;Length[{c}]>2

];

generateAction[a___]:=Message[generateAction::syntax,a];



(* create interaction list from action *)

(* if action is already given as list; can happen when using setSourcesZero manually *)
getInteractionList[action_List]:=action;

getInteractionList[action_]:=Module[
{},
Cases[action, S[__], \[Infinity]] /. S[a__] :> First /@ {a}	
]

getInteractionList[a___]:=Message[getInteractionList::syntax,a];



(* function to determine signs from fermions *)
(* Signs are represented by sf[a,b] where a is a single field and b a list of fields.
If a is a Grassmann field, a minus sign appears for each Grassmann field in b.
sf basically represent the sign introduced by moving the field a beyond the fields b.
*)

(* no fields left *)
sf[a_, {}] := 1
   
(* if the first field is not Grassmann, then no sf is needed *)
sf[a_, b___List] /; 
  bosonQ[a] || complexFieldQ[a] || antiComplexFieldQ[a] := 1

(* delete double fields *)
sf[a_, b___List] /; Sort@b =!= Union[b] := 
 sf[a, b /. {m1___, n_, m2___, n_, m3___} :> {m1, m2, m3}]
 
(* delete non-Grassmann fields *)
sf[a_, b___List] /; 
  Not[FreeQ[
    b[[All, 1]], _?bosonQ | _?complexFieldQ | _?antiComplexFieldQ]] :=
  sf[a, DeleteCases[b, _?bosonQ | _?complexFieldQ | _?antiComplexFieldQ]]
  
(* if elements of op are sorted, do not touch sf expressions *)
Sort[sf[a_, b_]] ^:= sf[a, b];

(* make signs explicit *)
getSigns[exp_] := 
 exp /. sf[f_?grassmannQ, l__] :> (-1)^Count[l, _?grassmannQ]




(* ::Section:: *)
(* Differentiation *)


(* for the first derivative differentiate once and then replace the fields by the corresponding expressions
using replaceFields *)

replaceFields[a__]:=replacementCalc@firstDerivReplacement[a]//Expand;


firstDerivReplacement[a_Times|a_Plus]:=firstDerivReplacement[#]&/@a;

firstDerivReplacement[a_?NumericQ]:=a;

firstDerivReplacement[op[b___,c_List]]:=op[Sequence@@(replacedField/@{b}),c];


(* no action on S, sf and CO *)
replacedField[a_S]:=a;
replacedField[a_CO]:=a;
replacedField[a_sf]:=a;
replacedField[a_?NumericQ]:=a;


replacementCalcStep[a_Times|a_Plus]:=replacementCalcStep/@a;

replacementCalcStep[a_?NumericQ]:=a;

replacementCalcStep[a_op]/;FreeQ[a,replacedField]:=a;

(* the utmost right two terms, quite simple *)
replacementCalcStep[op[a___,replacedField[{Q_,q_}],c_List]]:=(op[a,{Q,q},c]+
	op[a,sf[{Q,q},{{Q,q}}],P[{Q,q},c/.{d_,e_}:>{d,e}]]);

(* all higher terms *)
replacementCalcStep[op[a___,replacedField[{Q_,q_}],c__/;FreeQ[{c},replacedField]]]:=Module[{ind1,ind2},
	op[a,{Q,q},c]+
(* propagator and derivative w.r.t. field *)Plus@@((op[a,sf[{Q,q},{{Q,q}}], P[{Q,q},#],Sequence@@DeleteCases[{c},#1]])&)/@Cases[{c},{_?fieldQ,_}]+
(* propagator and derivative w.r.t. vertex *) Plus@@(op[a,P[{Q,q},{$dummyField,ind1=insDummy[]}],Sequence@@DeleteCases[{c},#1],derivVertex[#1,{$dummyField,ind1}]]&)/@{c}+
(* propagator and derivative w.r.t. propagator *) $signConvention Plus@@(op[a, P[{Q,q},{$dummyField,ind2=insDummy[]}],
	 Sequence@@DeleteCases[{c},#1],derivPropagator[#1,{$dummyField,ind2}]]&)/@{c}
];


replacementCalc[a_]:=FixedPoint[replacementCalcStep,a,50];


(* auxiliary function for plugging in the field into vertices and propagators *)
(* called with no fields to plugin *)
plugInFieldsV[{},b_V]:=b;
plugInFieldsV[{Q_?fieldQ,q_},b_V]:=Prepend[b,{Q,q}];
plugInFieldsV[newField_List,b__]:=plugInFieldsV[Most@newField, plugInFieldsV[newField[[-1]], b]];
 

(* when performing the differentiation a set of rules is applied in deriv;
this then calls derivAll which uses derivField, derivVertex and derivPropagators; 
the latter gets the correct factors if there are more of the same propagators using derivPropagator; *)


derivAll[a_,{Q_,q_}]/;Not@FreeQ[a,{Q,q},Infinity]:=Message[derivAll::index,{Q,q}];

derivAll[op[fvp___],{Q_,q_}]:=Module[{allShifts,i, permSign},
	
	
	allShifts=op[fvp];
	
	(* include signs from permuting the derivative through all expressions left of the target *)
	permSign[leftFvp_,{R_,r_}]:=sf[{R,r}, Cases[leftFvp, S[___] | P[__] | V[__] | {_?fieldQ,_}]/.{S:>Sequence, P:>Sequence, V:>Sequence}];
	
	derivField[op[fvp],{Q,q}]+
   + $signConvention Plus @@ Table[ReplacePart[op[fvp], i -> Sequence[permSign[Take[{fvp}, i-1], {Q,q}], derivPropagator[{fvp}[[i]], {Q, q}]]], {i, 1, Length[{fvp}]}]
   + Plus @@ Table[ReplacePart[op[fvp], i -> Sequence[permSign[Take[{fvp}, i-1], {Q,q}], derivVertex[{fvp}[[i]], {Q, q}]]], {i, 1, Length[{fvp}]}]
	
]	
	
	
(* if called with no derivatives given *)
deriv[a_]:=a;

deriv[a_,firstInd_,otherInds__]:=deriv[deriv[a,firstInd],otherInds];

deriv[a_,{Q_,q_}]:=sortDummies@Expand[a/.op[b__]:> derivAll[op[b],{Q,q}]];

deriv[a___]:=Message[deriv::syntax,a];


(* here the single op functions are orderd such that the derivatives can be done without sign problems *) 

derivRGE[a_,firstInd_,otherInds__]:=derivRGE[derivRGE[a, firstInd],otherInds];

derivRGE[a_,{Q_,q_}]:=Expand[a/.op[b__]:> derivAllRGE[op[b],{Q,q}]];


(*derivRGE[a___]:=Message[derivRGE::syntax,a];*)


derivAllRGE[op[fvp___], {Q_,q_}]:=Module[{i, permSign},

	(* TODO: Add deriv of fields *)
	(* include signs from permuting the derivative through all expressions left of the target *)
	permSign[leftFvp_,{R_,r_}]:=sf[{R,r}, Sequence@@@Cases[leftFvp, P[__] | V[__] | {_?fieldQ,_}]];
	
	$signConvention Plus @@ Table[ReplacePart[op[fvp], i -> Sequence[permSign[Take[{fvp}, i-1], {Q,q}], derivPropagator[{fvp}[[i]], {Q, q}]]], {i, 1, Length[{fvp}]}]
   + Plus @@ Table[ReplacePart[op[fvp], i -> Sequence[permSign[Take[{fvp}, i-1], {Q,q}],  derivVertex[{fvp}[[i]], {Q, q}]]], {i, 1, Length[{fvp}]}]
];

 

 
derivField[op[fvp___],{A_,i_}]:=Module[{
fields,
rest,
pos,
posToReplace},
	
	(* split into fields and the rest *)
	fields = Cases[{fvp},{_?fieldQ,_}];
	rest = Complement[{fvp},fields];
	(* get the positions of the corresponding fields *)
	pos=Position[fields,{A,_},1];
	(* factor for multiple fields; replace the first/last (anti-fermions/fermions) field by the one w.r.t. which we differentiate and delete the same in the external fields; *)
	If[pos!= {} (* check if there are fields *),
		(posToReplace:=pos[[1]];
		
		(* add sf for anticommutations with fermions *)
		Length@pos op[sf[{A,i},Take[fields,posToReplace[[1]]-1]], Sequence@@(rest/.fields[[posToReplace[[1]]]]:> {A,i}),Sequence@@Delete[fields,posToReplace]]),
		0 (* no fields *)
		]
];

derivField[op[S_],{Q_,q_}]:=0;



derivVertex[V[fields__],{Q_,q_}]:=V[{Q,q},fields];

derivVertex[S[fields__],{Q_,q_}]:=0;

derivVertex[P[a__],{Q_,q_}]:=0;

derivVertex[CO[a__],{Q_,q_}]:=0;

derivVertex[a_List,{Q_,q_}]:=0;

derivVertex[a_?NumericQ,{Q_,q_}]:=0;

derivVertex[a_dR,{Q_,q_}]:=0;

derivVertex[a_sf,{Q_,q_}]:=0;
 

derivPropagator[V[a__],{Q_,q_}]:=0;

derivPropagator[CO[a__],{Q_,q_}]:=0;

derivPropagator[a_List,{Q_,q_}]:=0;

derivPropagator[a_?NumericQ,{Q_,q_}]:=0;

derivPropagator[S[__],{Q_,q_}]:=0;

derivPropagator[sf[_,_],{Q_,q_}]:=0;

derivPropagator[P[field1_List,field2_List],{Q_,q_}]:=ReleaseHold@Module[{dummy1,dummy2},
	Hold@Sequence[sf[{Q,q},{field1,{$dummyField,dummy1=insDummy[]}}],P[field1,{$dummyField,dummy1}],V[{Q,q},{$dummyField,dummy1},{$dummyField,dummy2=insDummy[]}],
		P[{$dummyField,dummy2},field2]]
];


derivPropagators[a_op,{Q_,q_}]:=Module[
{props,propConnections,posInd,represConnections,connectionList,propFactors,propagatorId},

(* test if two propagators can be considered equal *)
propagatorId[b_,c_]:=(b[[2]]===c[[2]](* connections have to be the same *)&&Cases[Sort@b[[1]],{_?fieldQ,_}][[All,1]]===Cases[Sort@c[[1]],{_?fieldQ,_}][[All,1]](* fields have to be the same *));

(* get all propagators *)
props=Cases[a,_P];

(* get vertices they connect *)
propConnections=Function[{prop},
posInd=Position[a,Alternatives@@prop];
{prop,DeleteCases[a[[Sequence@@#[[1]]]]&/@posInd,_P]}]/@props;

(* delete appearances of sf function *)
propConnections = propConnections/.sf[_,_]:>Sequence[];

(* create a list of the factor, 1 representative and the connected vertices *)
represConnections=Union[propConnections,SameTest->propagatorId];

propFactors=Function[prop,Length@Select[propConnections,propagatorId[#,prop]&]]/@represConnections;

connectionList=Transpose[{propFactors,represConnections}];

ReleaseHold[Apply[Plus,(a/.#[[2,1]]:> #[[1]]Hold@derivPropagator[#[[2,1]],{Q,q}])&/@connectionList]]
];


(* insert a regulator; don't worry about the fields, they are fixed since the sources should already be zero *)
(* a minus sign appears here (no $signConvention because no vertices are involved) *)
insertRegulator[P[{field1_,ind1_},{field2_,ind2_}]]:=ReleaseHold@Module[{dummy1,dummy2},
	Hold@Sequence[-dR[{field2,dummy1=insDummy[]},{field1,dummy2=insDummy[]}],
		P[{field1,ind1},{field2,dummy1}],P[{field1,dummy2},{field2,ind2}]]
];

(* differentiate with respect to t tilde; this amounts to replace a propagator D by D dt R D *)
derivPropagatorsdt[a_op]:=Module[
{props,propConnections,posInd,represConnections,connectionList,propFactors,propagatorId},

(* test if two propagators can be considered equal *)
propagatorId[b_,c_]:=(b[[2]]===c[[2]](* connections have to be the same *)&&Cases[Sort@b[[1]],{_?fieldQ,_}][[All,1]]===Cases[Sort@c[[1]],{_?fieldQ,_}][[All,1]](* fields have to be the same *));

(* get all propagators *)
props=a[[Sequence@@#]]&/@Position[a,_P];

(* get vertices they connect *)
propConnections=Function[{prop},
posInd=Position[a,Alternatives@@prop];
{prop,DeleteCases[a[[Sequence@@#[[1]]]]&/@posInd,_P]}]/@props;

(* create a list of the factor, 1 representative and the connected vertices *)
represConnections=Union[propConnections,SameTest->propagatorId];

propFactors=Count[propConnections[[All,2]],#]&/@represConnections[[All,2]];
propFactors=Function[prop,Length@Select[propConnections,propagatorId[#,prop]&]]/@represConnections;

connectionList=Transpose[{propFactors,represConnections}];

ReleaseHold[Apply[Plus,(a/.#[[2,1]]:> #[[1]]Hold@insertRegulator[#[[2,1]]])&/@connectionList]]
];




(* ::Section:: *)
(* Identification *)


(* for RGEs there is an algorithm that will always work, since RGEs are one-loop only;
note that the coefficient can become 0 if the signs are wrong, e.g., -1/2+1/2=0 *)

(* in case there is only one term *)
identifyGraphs[exp_Times,extFields_List]:=exp;

identifyGraphs[exp_op,extFields_List]:=exp;

identifyGraphs[a_?NumericQ,opts___]:=a;


(* more terms *)
identifyGraphs[exp_Plus,extFields_List]:=Module[{classes,ops,equalOps,orderedExp},

orderedExp=sortCanonical[exp, extFields];

(* split off the numerical factors; syntax: {{factor1,op1},{factor2, op2},{factor3,op3},...} *)
ops=Replace[List@@Expand@orderedExp/.Times[b_?NumericQ,c_]:> {b,c},d_op:> {1,d},{1}];

(* identify only within single classes to improve performance;
classes are identified by 
-) number and types of propagators
-) types of vertices and their external legs *)
classes=Flatten[GatherBy[ops, {
  Sort@Cases[#1, P[q1_, q2_] :> {q1[[1]], q2[[1]]}, 2]&,
  Sort[(Cases[#, V[__]|CO[__], 2]/. {Q_?fieldQ, q_} :> {Q} /; 
      Not@MemberQ[extFields[[All, 2]], q])/.V[a__]:> Sort[V[a]]/.CO[a__]:>Sort[CO[a]]
      ] &}],1]; 
      
(* add up equivalent graphs *)
equalOps=Flatten[#//.{b___,c_List,d___,e_List,f___}:> {b,{c[[1]]+e[[1]],c[[2]]},d,f}
	/;getGraphCharacteristic[c[[2]],extFields]===getGraphCharacteristic[e[[2]],extFields]&/@classes,1];

(* multiply with numerical factors *)
equalOps[[All,1]].equalOps[[All,2]]

];

identifyGraphs[a___]:=Message[identifyGraphs::syntax,a];


(* get  characteristic of a graph
A characteristic of a graph is a tree-like structure containing information about all the neighbouring vertices.
The first entry is the vertex with the first external field.
The second entry are its neighbours:
The connecting field and the vertex, but again in list form, so that they contain information about their neighbours recursively.
Further vertices only contain the new neighbour.
At the end all indices except the external ones are removed and the vertex arguments are sorted.
All equal graphs have then the same characteristic. *)

(* auxiliary function to determine the neighbours *)
(* if there are no external fields, it must be a zero-point RGE --> no starting vertex, just do canonial sorting  *)
getNeighbours[exp_, {}] := sortCanonical[exp, {}];

(* first instance: start at the vertex with the first external index *)
getNeighbours[exp_, allExtFields_List] := 
  Module[{mainVertex, extFieldsOfMainVertex, connectingFields, allExtFields2, fieldValues, i, sortByExtField},
  	
  (* add the real external fields to the external legs (given by allExtFields); in this way the connecting fields are determined correctly *)
   allExtFields2=Cases[exp,{Q_,q_}/;Count[exp,{Q,q},\[Infinity]]==1,{2}];

   (* determine the first vertex *)
   mainVertex = Select[exp, Not@FreeQ[#, allExtFields[[1]]] &][[1]];
   (* if there is a single loop that starts and ends at the starting vertex, we have to delete that *)
   mainVertex = mainVertex/. V[a___, b_, c___, b_, d___] :> V[a, b, c, d]; 
   
   (* the external fields in the starting vertex; could be others too besides the first external field *)
   extFieldsOfMainVertex = 
    Cases[mainVertex, Alternatives @@ (# & /@ allExtFields2)];

   (* determine the internal fields *)
   connectingFields = 
    List @@ Complement[mainVertex, 
      Sequence @@ (V[#] & /@ extFieldsOfMainVertex)];
   connectingFields = List@@mainVertex /. {f_?fieldQ, j_}:>Sequence[]/;MemberQ[allExtFields, {f,j}];
   (* result: {starting vertex, {{field which connects to neighbour 1, neighbour 1},{field which connects to neighbour 2, neighbour 2}}} *)
   
   fieldValues = Association @@ Table[allExtFields[[i, 2]] -> i, {i, 1, Length[allExtFields]}];
   
   (* sort the connections by the external fields *)
	   (*sortByExtField[a_]:=SortBy[
	      a, Function[b,Select[b, Not@FreeQ[#, Alternatives @@(allExtFields[[All,2]])]&]]];*)
   sortByExtField[a_]:=a;(* use sortCanonical before, so the vertex arguments should already be sorted *)
   
   (* TODO: commented the second FreeQ to work with diagrams with self-loops, 
      	unclear if that works for other cases, because it was introduced on purpose; up to now no problems *)
    {mainVertex, 
     Function[cF,
      	{cF, Select[exp, Not@FreeQ[#, cF] (*&& FreeQ[#, allExtFields[[1]]]*) &][[1]]}] /@ connectingFields}
   ];
   

(* neighbour for further instances *)
getNeighbours[exp_, startingField_List, v_V, allExtFields_List] := Module[{connectingFields, fieldsDone},

(* determine fields which yield no new neighbours *)
   fieldsDone = 
    Cases[List @@ v, 
     Alternatives @@ Union[allExtFields, {startingField}]];
  
(* get new neighbour, delete fields already done *)
   (*connectingFields = List @@ Complement[v, V @@ fieldsDone];*)
   connectingFields = List @@ (v /. (Rule[#, Sequence[]] &/@fieldsDone));
   
(* result: {vertex, {field which connects to new neighbour, new neighbour}} *)
   {v, Function[
      cF, {cF, 
       (Select[exp, 
         Not@FreeQ[#, cF] && FreeQ[#, Alternatives@@fieldsDone] &]/.{}:>{{}}(* avoids problems with external fields, which have no connections *))[[1]]}] /@ 
     connectingFields} 
];


(* special case: tadpole like diagram -> only one vertex *)
getGraphCharacteristic[graph_op, extFields_List] /; 
   Count[graph, V[__]|CO[___]|S[___]] == 1 :=(* getGraphCharacteristic[graph, extFields] =*)
 graph /. P[__] :> Sequence[] /. {Q_Symbol, q_Symbol} :> {Q} /; 
      Not@MemberQ[extFields[[All, 2]], q] /. V[a__] :> Sort[V[a]];

getGraphCharacteristic[graph_op, extLegs_List] := 
 Module[{props, verts, id,firstNeighbours, allFieldsInV, extFields, extFieldsRotated},

  props = Cases[graph, P[__]];
  verts = Cases[graph, V[__]|S[___]|CO[___]]/.S[a___]:>V[a]/.CO[a___]:>V[a];(* rewrite to vertices V for getNeighbours *)
  
  (* transform into vertices only; these rules produce apparently non-
  existent vertices, but the connections can still be identified, 
  which is the important part *);  
  id = verts /. (props /. P :> Rule);
  
  (* at this point we can determine all legs connected to an external field or a derivative; before this was not possible, because the external fields appeared twice (in V and op) *)
  allFieldsInV=Cases[id,{_?fieldQ,_},\[Infinity]];
  extFields=Select[allFieldsInV,Count[allFieldsInV,#]==1&];
  (* rotate extFields until the first derivative is on position 1; this is necessary so that equal graphs with opposite directions of the fields can be identified  *)
  (*extFieldsRotated=FixedPoint[RotateLeft[#]&, extFields, Length@extFields, SameTest->(Not[FreeQ[#2[[1]], extLegs[[1]]]] &)];*)
  extFieldsRotated=extLegs;
  (* find neighbours from there *)
  firstNeighbours = getNeighbours[id, extFieldsRotated];
 
  (* repeat the process until the loop is closed *)
  Nest[(# /. {a_, b_V|b_CO} :> {a, getNeighbours[id, a, b, extFieldsRotated]} )&,
     firstNeighbours,	
     Max[0, Floor[(Length@props-2)/2]]]
    /. {Q_?fieldQ, q_} :> {Q} /; 
      Not@MemberQ[extFieldsRotated[[All, 2]], q] /. V[a__] :> Sort[V[a]]
    (*/. {a_, {b_List, c_List}} :> {a, Sort[{b, c}]}*) (* sorting the first neighbours;
    	this identifies graphs with the opposite ordering direction of the external legs *)
]


(* Sorts expressions into a canonical order of field arguments. 
Includes effect of the deprecated function orderFermions.
External fields are not fully included here.*)

(* if there are dummy fields in the expression, sorting does not work *)
sortCanonical[b_op, derivatives_List]/;Not[FreeQ[b,$dummyField]]:=(Message[sortCanonical::dummyFields, b];b)

sortCanonical[b_op, derivatives_List] := 
 Module[{ordered, fieldValues, i, intIndices, extFieldIndices, props, const,
 	connectedLeg, intIndexAss, intVertsAss, orderV, orderP, orderR, intVerts, extVerts,
 	complexPropRules,extFields},

  (* constant to distinguish internal from external indices *)
  const = 10^10;
  
  (* rules to sort complex propagators *)
  complexPropRules = {P[{phi_?fieldQ, ii1_}, {phibar_?fieldQ, ii2_}] /; (complexFieldQ[phibar] && antiField[phi]===phibar) :> P[{phibar, ii2}, {phi, ii1}]};
  
  (* extract external fields *)
  extFields = Cases[b, {_?fieldQ,_}];
  
  (* extract propagators *)
  props =  Cases[b, P[___]];
  
  (* get vertices with external legs *)
  extVerts = Cases[b, V[a__]|S[a__]|CO[a__]/;Not@FreeQ[{a}, Alternatives@@derivatives[[All,2]]] && Length[{a}]>2];
  
  (* get vertices without external legs *)
  intVerts = Cases[b, V[a__]/;FreeQ[{a}, Alternatives@@derivatives[[All,2]]]];
  
  (* assign each index a number for sorting: 
  external ones by position in derivatives, 
  internal ones by their connection to external ones*)

  extFieldIndices = Cases[b, {a_?fieldQ, j_}:>j];
  
  (* get internal indices *)
  intIndices = Complement[Union@Cases[b, {a_?fieldQ, j_} :> j, Infinity], derivatives[[All, 2]]];
  
  (* delete external field indices *)
  intIndices = Complement[intIndices, extFieldIndices];
  
  (* set up association for external indices *)
  fieldValues = Association @@ Table[derivatives[[i, 2]] -> i, {i, 1, Length[derivatives]}];
  
  (* determine index that connects the index fieldInd with a vertex *)
  connectedLeg[fieldInd_]/;MemberQ[extFields[[All,2]],fieldInd]:={};
  connectedLeg[fieldInd_] := DeleteCases[Select[List @@@ props, Not@FreeQ[#, fieldInd] &][[1]], {_, fieldInd}][[1, 2]];
  
  (* set up association for internal indices: 
  extract vertices that internal indices connect to, 
  then determine their association values by the smallest value of an external field there and add a constant *)
   
  (* assign internal vertices a value for ordering *)
  (* extract for each internal index the vertex it connects to *)
  intIndexAss = {#, (Cases[b, V[a__]|S[a__]|CO[a__] /; Not@FreeQ[{a}, connectedLeg@#]]/.{}:>{{}})[[1]]} & /@ intIndices;
  (* remove internal vertices *)
  intIndexAss = Select[intIndexAss, Not@MemberQ[intVerts, #[[2]]]&];
  (* assign a field value based on the lowest field value of the external legs of the connected vertex *)
  intIndexAss = intIndexAss /. V[a___]|S[a___]|CO[a__] :> const + Sort[(fieldValues[#[[2]]] & /@ Select[{a}, MemberQ[derivatives, #] &])][[1]];
  
  (* combine associations *)
  fieldValues = Join[fieldValues, Association @@ Rule @@@ intIndexAss];
  
  (* handle internal vertices by assigning a value based on the connected vertices *)
  (* take their indices *)
  intVertsAss = (List@@#)[[All,2]] & /@ intVerts;
  (* get all vertices they connect to and take their indices only*)  
  intVertsAss = Map[connectedLeg, intVertsAss, {2}];
  intVertsAss = Map[Function[ind, Select[extVerts, Not@FreeQ[#, ind] &]], intVertsAss, {2}] /. {a_?fieldQ, 
   j_} :> j /. S:>List /. V:>List /. CO:>List;
  (* get field values for all indices *)
  intVertsAss = Map[fieldValues, intVertsAss, {4}];
  (* determine lowest field value of vertex that internal indices connect to (those with missing keys) *)
  intVertsAss = Map[Cases[#, Missing["KeyAbsent", ind_] :> {ind, Min@Select[#, NumericQ]}] &, intVertsAss, {3}];
  (* assign field value *)
  intVertsAss = Union@Flatten[Map[Rule[#[[1]], 2*const + #[[2]]] &, intVertsAss, {4}]];
  (* do the same for the connecting counterparts *)
  intVertsAss = Flatten[Join[intVertsAss, intVertsAss/.(j_->ind_):>(connectedLeg[j]->ind)]]; 
   
  (* add the new associations to the field value function *)
  fieldValues = Join[fieldValues, Association @@ intVertsAss];
    
  (* order fields of a vertex as bosons, antiFermions, 
  fermions and internally by fieldValues *)
  
  orderV[V[a__]] := V[Sequence @@ Join[
  	SortBy[Select[{a}, cFieldQ[#[[1]]] &], fieldValues[#[[2]]] &],
  	SortBy[Select[{a}, antiFermionQ[#[[1]]] &], fieldValues[#[[2]]] &], 
    SortBy[Select[{a}, fermionQ[#[[1]]] &], (fieldValues[#[[2]]]) &]]];
  orderV[S[a__]] /; Length[{a}]>2 := orderV[V[a]]/.V:>S;
  orderV[CO[a__]] /; Length[{a}]>2 := orderV[V[a]]/.V:>CO;
  orderV[S[a__]] := S[a];
  orderV[CO[a__]] := CO[a];
  
  orderP[P[{Q1_?antiFermionQ, q1_}, {Q2_?fermionQ, q2_}]] := P[{Q2, q2}, {Q1, q1}];
  orderP[P[{Q1_, q1_}, {Q2_, q2_}]] := P[{Q1, q1}, {Q2, q2}];
  
  orderR[dR[{Q1_?fermionQ, q1_}, {Q2_?antiFermionQ, q2_}]] := dR[{Q2, q2}, {Q1, q1}];
  orderR[dR[{Q1_, q1_}, {Q2_, q2_}]] := dR[{Q1, q1}, {Q2, q2}];
    
  ordered = b /. V[a__] :> orderV[V[a]] /. S[a__] :> orderV[S[a]] /. CO[a__] :> orderV[CO[a]] 
  	/. P[a__] :> orderP[P[a]] /. dR[a__] :> orderR[dR[a]];
  
  (* get signature sign of original and ordered expression for the relative sign *)
  getSignature@b getSignature@ordered ordered/.complexPropRules
]


(* apply to a complete expression *)
sortCanonical[b_, derivatives_List] := b/. op[a___]:> sortCanonical[op[a], derivatives]


(* Auxiliary function to get a sign for an expression. Used to compare with reordered expression to get the relative sign. The obtained sign itself has no meaning. *)
getSignature[a_op] := Module[{verts},
  (* get vertices and keep only Grassmann fields *)
  verts = DeleteCases[
    Cases[a, V[___]|S[___]|CO[___]|P[___]|dR[___]], {_?bosonQ, _} | {_?complexFieldQ, _} | {_?antiComplexFieldQ, _}, {2}];
  (* return sign *)
  Times @@ (Signature /@ verts)
]




(* ::Section:: *)
(* Set Sources to Zero *)

(* standard test checking for conservation of Grassmann numbers
here we assume that only Grassmann numbers of fields are conserved that have a propagator *)
grassmannTest[a_V,fields_List] := Module[{fermionsList},

fermionsList=Cases[fields,{_?(fermionQ@#||antiFermionQ@#&),_?(fermionQ@#||antiFermionQ@#&)}];

And@@(Count[a,#[[1]],\[Infinity]]===Count[a,#[[2]],\[Infinity]]&/@fermionsList)
];


(* standard test checking if fields that have only even interactions only appear in vertices with even legs of this field;
note that complex fields are tested in grassmannTest *)

evenBosonTest[a_V,evenBosons_List,evenComplexBosons_List] := Module[{},

(* test if complex bosonic fields and even (real) bosonic fields are even in the vertices *)
And@@(Count[a,#[[1]],\[Infinity]]===Count[a,#[[2]],\[Infinity]]&/@evenComplexBosons) &&
And@@(EvenQ@Count[a,#,\[Infinity]]&/@evenBosons)

];


(* RGE derivations use a different convention for Grassmann fields; this is realized via special propagator replacement rules and differentiation rules of propagators;
the former is decided here by using setSourcesZeroRGE, which just sets the option propagatorCreationRules to RGERules but otherwise uses the normal setSourcesZero;
createPropagatorRules then knows what to do;
note that the RGE convention is considered superior as it is unique, while for the DSE convention there are still problems for vertices and two loops, like, e.g., the A c cb vertex
-> it is recommended to order Grassmann and anti-Grassmann fields starting with the second derivative, so use A cb c here;
unfortunately not known/possible how to use the RGE convention for DSEs *)
setSourcesZeroRGE[all___,opts___?OptionQ]:=
	setSourcesZero[all,propagatorCreationRules->RGERules,opts];



(* the following two calls of setSourcesZero replace the vertices appropriately;
this is done here in order to exclude double counting, because the real code of setSourcesZero is called for every field *)
(* interactions given *)
setSourcesZero[a_,L_,extLegs_List,ownAllowedPropagators_List:{},vertexTest_Symbol|vertexTest_Function,opts___?OptionQ]:=Module[
	{fields,numberExtFields,verticesAtMin,interactions,truncation},
	
	(* get interactions from the action *)
	interactions=getInteractionList@L;

	(* list of different fields *)
	fields=Union@Flatten@interactions;
	
	(* for DSEs in the broken phase we need an ansatz for the effective action; in order to allow abbreviations like {phi,6} create a generated action and extract the interactions from it *)
	truncation=getInteractionList@generateAction[(ansatz/.Join[{opts},Options@doDSE])];

	(* number of external fields allowed;
	the maximal number for a vertex is given by the highest interaction in the action -3 (for 1 ext. and 2 internal legs);
	for RGEs we have the ansatz given in the action;
	for DSEs the ansatz is supplied with the option ansatz *)    	
    numberExtFields=If[(symmetry/.Join[{opts},Options@doRGE])===broken,
    	If[(propagatorCreationRules/.{opts})=!=RGERules,Max[Max[Length/@truncation]-3,0(* make sure it is not negative *)],
    	   Max[Max[Length/@interactions]-3,0(* make sure it is not negative *)],0],
    		0,
    		0];
    		
    (* replace the vertices; truncations are done later *)
    verticesAtMin=replaceVerticesBrokenPhase[a,fields,numberExtFields,interactions];
 
    (* set sources to zero field by field *)
	setSourcesZeroDo[verticesAtMin,L,extLegs,truncation,ownAllowedPropagators,numberExtFields,vertexTest,opts]
];

(* if no vertex test given *)
setSourcesZero[a_,L_,extLegs_List,ownAllowedPropagators_List:{},opts___?OptionQ]:=setSourcesZero[a, L, extLegs, ownAllowedPropagators, (True&), opts]

setSourcesZero[a___]:=Message[setSourcesZero::syntax,a];



setSourcesZeroDo[a_,L_,extLegs_List,truncation_List,ownAllowedPropagators_List:{},numberExtFields_Integer,vertexTest_Symbol|vertexTest_Function,opts___?OptionQ]:=Module[
{test,standardTest,propRules,propReplaced,allowedPropagators,fields,evenBosons,evenComplexBosons,bosons,
verticesAtMinTrunc, verticesAtMinTruncMarked,vertexPatterns,vertexTestTruncation,interactions,truncationTest},


(* for DSEs in the broken phase this test is required; it compares a vertex against the truncation ansatz given; can also be used in the symmetric phase *)
truncationTest[v_V, {}]:=True;
truncationTest[v_V, ans_]:=MemberQ[Sort/@ans,List@@Sort[v[[All,1]]]];
	
(* get interactions from the action *)
interactions=getInteractionList@L;

fields=Union@Flatten@interactions;

(* get allowed propagators from the action or use the user-defined allowed propagators *)
(*allowedPropagators=Join[ownAllowedPropagators,Cases[interactions, b_List /; Length@b == 2]];*)
allowedPropagators = ownAllowedPropagators/.{}:>Cases[interactions, b_List /; Length@b == 2];

(* if there are fields that have only even interactions there cannot be a vertex with an odd number of these fields;
for fermions this check is done with grassmannTest *)
bosons=Select[Union@Flatten@interactions, cFieldQ];
(* add the fields which can be specified by the user as even, when one uses doRGE and no truncations *)
evenBosons=Union[Cases[userEvenFields/.Join[{opts},Options@doRGE],_?fieldQ],Function[boson, And @@ ((EvenQ@Count[#, boson] & /@ interactions)/.{}:>False)/. {True :> boson, False :> Sequence[]}] /@ bosons];
evenComplexBosons=Union[Cases[userEvenFields/.Join[{opts},Options@doRGE],{_?fieldQ,_?fieldQ}], Function[complexBoson, And @@ ((Count[#, complexBoson[[1]]]==Count[#, complexBoson[[2]]] & /@ interactions)/.{}:>False)/. {True :> complexBoson, False :> Sequence[]}]/@({#,antiField[#]}&/@Select[bosons, complexFieldQ])];

(* determine the possible interactions from the truncation given as input;
use the DSE function generateAction to create bare vertices and use them as patterns in a test;
only applies for RGEs *)
vertexPatterns=Cases[L,S[b1_,b2_,b3__]:>V[b1,b2,b3](* ensure it's a vertex and not a two-point fct. *),\[Infinity]] /. {Q_?fieldQ, _} -> {Q, _};
vertexTestTruncation[v_V,{},___]:=True; (* no truncation given, so take everything *)
vertexTestTruncation[v_V,vPs_,RGERules]:=MatchQ[Sort@v,Alternatives@@(Sort/@vPs)]; (* take all vertices of truncation if we have an RGE; sorting is required since the fields come in no special order *)
vertexTestTruncation[v_V,___]:=True;

(* define a test that is used when vertexTest is empty *)
standardTest[b_]:=evenBosonTest[b,evenBosons,evenComplexBosons];

(* test function for vertices *)
test[b_]:=And[
	(* grassmannTest is only required for DSEs, since the ansatz for the action fully determines the allowed vertices for RGEs;
	thus switched off for RGEs to avoid problems with mixed fermions *)
	If[(doGrassmannTest/.Join[{opts},Options@setSourcesZero]) &&
		(propagatorCreationRules/.Join[{opts},Options@setSourcesZero])===DSERules,
		grassmannTest[b,allowedPropagators],True,True], (* conserve Grassmann number *)
	(* truncationTest now used always *)(*If[(propagatorCreationRules/.Join[{opts},Options@doDSE])=!=RGERules&&(symmetry/.Join[{opts},Options@doDSE])===broken,truncationTest[b,truncation],True,True], (* DSEs, broken phase, compare against truncation ansatz *)*)
	truncationTest[b,truncation],
	If[(propagatorCreationRules/.Join[{opts},Options@setSourcesZero])===RGERules,truncationTest[b,interactions],True,True], (* RGEs, always compare against interactions in the action *)
	standardTest[b], (* conserve boson number if they only have even interactions *)
	vertexTest[b](* user-defined test *)
];

propRules=createPropagatorRules[allowedPropagators,fields,opts];

(* truncation: only allow certain number of external fields;
in DSE calculations it should be possible to set only a certain type of fields to zero, so do not truncate if we are in the symmtric phase;
numberExtFields gives the external fields per vertex;
add the marker ext to all external fields, then count them per vertex and set vertices with too many fields to zero, finally remove the marker *)
verticesAtMinTruncMarked=a/.op[e___,{F_?fieldQ,f_},g___]:>(op[e,{F,f},g]/.{F,f}:>{{F,f},ext});
verticesAtMinTrunc=verticesAtMinTruncMarked/.exp_V|exp_S|exp_CO:>0/;Count[exp,{_List,ext}]>numberExtFields/.{l_List,ext}:>l;

(* apply replacement rules for the propagators; note that here we can get a sum of ops;
delete propagators that are not allowed *)
propReplaced=verticesAtMinTrunc/.propRules/.P[{q1_,_},{q2_,_}]:>0 /; Not@MemberQ[allowedPropagators, {q1, q2}|{q2,q1}];

(* for each term in sum replace the vertices separately; 
close the indices *)
Expand@replaceVertices[propReplaced,test,extLegs,opts]

];



replaceVertices[c_Times,rest___]:=replaceVertices[#,rest]&/@c;

replaceVertices[c_Plus,rest___]:=replaceVertices[#,rest]&/@c;

replaceVertices[a_?NumericQ,___]:=a;

replaceVertices[d_op,vertexTest_Symbol|vertexTest_Function,extFields_List,opts___?OptionQ]:=Module[{propagators,propInds, propagatorsF,
	propagatorsB, propIndsF, propIndsB,c,vertsReplaced},

(* indices of the propagators *)
c=d/.traceIndex1:>traceIndex2 (* close the trace here so that all variants are taken into account *);
propagators=Cases[c,P[___],Infinity];
propInds=Cases[propagators,{_?fieldQ,_},Infinity];

propagatorsF=Select[propagators,(Head@#[[1,1]]==fermion)&];
propagatorsB=Select[propagators,Head@#[[1,1]]==boson&];

propIndsF=Cases[propagatorsF,{_?fieldQ,_},Infinity];
propIndsB=Cases[propagatorsB,{_?fieldQ,_},Infinity];

(* replace the indices of the vertices, delete the vertices that do not exist *)
vertsReplaced=c/.Plus:>List/.{{$dummyField,ind_}:> (Flatten@Cases[propInds,{_,ind}])}
	/.V[___,{},___]:>0/.V[b__]:> 0/;(Not@vertexTest[V[b]])
	
];




(* determine the vertex in the case that the field expectation value is not zero -> higher vertices multiplied by fields remain *)
(* give the vertex with external fields *)
newVertex[oldVertex_V, newFields_List] := 
 Module[{fields, dInds, indsFieldsCombined,symmetryFactor},
  
  (* lists of fields, the corresponding symmetry factors and dummy indices *)
  fields= Union@newFields;
  symmetryFactor=Times @@ Factorial /@ (Count[newFields, #] & /@ fields);
  dInds = Table[insDummy[], {Length@newFields}];
  
  (* combine fields and dummy indices *)
  indsFieldsCombined = Partition[Riffle[newFields, dInds], 2];
  
  (* plug fields into vertices and as external fields *)
  (* the order of the external fields is fixed later and of no meaning here *)
  ReleaseHold@{Sequence @@ indsFieldsCombined, 
    plugInFieldsV[indsFieldsCombined,oldVertex]/symmetryFactor}
];


(* sum up all possibilities until a given order *)
newVertexSum[exp_, allNewFields_List] := 
(exp /. {V[a___] :> newVertex[V[a], #]} & /@ allNewFields);
  

(* insert the new vertices; this complicated way is necessary because every new vertex requires its own dummy indices *)
(* in case the external fields are zero *)
replaceVerticesBrokenPhase[exp_,interactions_,0,interactions_List]:=(exp);

replaceVerticesBrokenPhase[exp_,fields_,numberExtFields_Integer,interactions_List]:=(*improves speed but may be unstable for too many calculations; replaceVerticesBrokenPhase[exp,fields,numberExtFields]=*)Module[
	{fieldsToAdd},

	(* TODO: this is a quick way around and calculates definitely too many vertices; still the speed is acceptable but not great -> should be improved by restricting to required vertices only*)
	
	(* create list of all field combinations and bring fermions and antifermions into the proper order;
	due to the way this list is created the order of different Grassmann fields is already correct, e,g. {antiF1, antiF2, F2, F1} *)
	fieldsToAdd = Flatten[Table[Union[Sort /@ Tuples[fields, i]], {i, 0, numberExtFields}], 1]//. {a___, b_, c___, d_, e___} :> {a, d, c, b, e} /; 
    fermionQ@b && antiFermionQ@d;
	
	(exp/.V[a___]:>(Plus@@op@@@newVertexSum[V[a],fieldsToAdd]))
]
 


(* for Grassmann fields ordering might be necessary to have anti-fields left of fields;
this can yield the minus sign of closed fermion loops *)
(* orderFermions is superseded by sortCanonical *)

orderFermions[a_Plus|a_Times]:=orderFermions/@a;

orderFermions[a_?NumericQ]:=a;

orderFermions[a_op] := Module[{fermions, bosons, antiFermions,orderF, orderB, orderExtFields, extFields, swapInnerFermions},

Message[orderFermions::superseded];

(* get lists of fermions, antiFermions and bosons from INTERNAL fields *)
fermions=Union@Cases[a,{f_?fermionQ,_}:>f,{2}];
antiFermions=Union@Cases[a,{f_?antiFermionQ,_}:>f,{2}];
bosons=Union@Cases[a,{f_?cFieldQ,_}:>f,{2}];

(* get external fields *)
extFields = Union@Select[Union@Cases[a, {_?fieldQ, _}, Infinity], Count[a, #, Infinity] == 1 &];

(* all bosons to the right and sort alphabetically to have a canonical order when using Feynman rules in doAE *)
orderB[b_op]:=b//.V[c1___,d1_,e1_,f1___]:>V[c1,e1,d1,f1]/;MemberQ[bosons,e1[[1]]]&&(MemberQ[fermions,d1[[1]]]||MemberQ[antiFermions,d1[[1]]])
	//. V[c2___, d2_, e2_, f2___] :> V[c2, e2, d2, f2] /; (Not@OrderedQ[{d2, e2}]&&MemberQ[bosons,e2[[1]]]&&MemberQ[bosons,d2[[1]]]);
 
(* all fermions to the right of the anti-fermions, order fermion names canonically, order by index *)
orderF[Times[c_,b_op]]:=c orderF[b];
orderF[b_op]/;fermions==antiFermions=={}:=b;
orderF[b_op]:= b //. 
  {V[e___, f_, g_, h___] :> -V[e, g, f, h] /;
    MemberQ[fermions, f[[1]]] && MemberQ[antiFermions, g[[1]]],
   (* S[a___, b_, c_, d___] :> -S[a, c, b, d] /; 
    MemberQ[fermions, c[[1]]] && MemberQ[antiFermions, b[[1]]],*)
    P[{Q1_?antiFermionQ,q1_},{Q2_?fermionQ,q2_}]:>-P[{Q2,q2},{Q1,q1}],
    dR[{Q1_?fermionQ,q1_},{Q2_?antiFermionQ,q2_}]:>-dR[{Q2,q2},{Q1,q1}]}/.
    V[e___?(Not@fermionQ[#[[1]]] &), f__?(fermionQ[#[[1]]] &), g___?(Not@fermionQ[#[[1]]] &)] :> Signature[{f}] V[e, Sequence@@(Sort[{f}]), g]/.
    V[e___?(Not@antiFermionQ[#[[1]]] &), f__?(antiFermionQ[#[[1]]] &), g___?(Not@antiFermionQ[#[[1]]] &)] :> Signature[{f}] V[e, Sequence@@(Sort[{f}]), g];
    
orderExtFields[Times[c_,b_op]]:=c orderExtFields[b];
orderExtFields[b_op]:=b
	(* sort all external fields to the left *)
	/. op[a1___, b1_List, c1_V | c1_P | c1_dR | c1_S, e1___] :> op[a1, c1, b1, e1]
	/. op[a2___, b2_List, c2___, d2_List, e2___] :> -op[a2, d2, c2, b2, e2] /; fermionQ@b2[[1]] && antiFermionQ@d2[[1]];

(* Change inner (only) fermion legs if they are not in canonical order. Treat fermions and antiFermions separately since that ordering should not be touched. *)
swapInnerFermions[Times[c_,b_op]]:=c swapInnerFermions[b];	
swapInnerFermions[exp_op] := 
 exp /. V[aa___, {b_?fermionQ, ib_}, {c_?fermionQ, ic_}, 
    d___] :> -V[aa, {c, ic}, {b, ib}, d] /; 
    Not@MemberQ[extFields,{b,ib}] && Not@MemberQ[extFields,{c,ic}] && Not@innerFermionsCanonicalQ[exp, V[aa, {b, ib}, {c, ic}, d]]
    /.V[aa___, {b_?antiFermionQ, ib_}, {c_?antiFermionQ, ic_}, 
    d___] :> -V[aa, {c, ic}, {b, ib}, d] /; 
    Not@MemberQ[extFields,{b,ib}] && Not@MemberQ[extFields,{c,ic}] && Not@innerFermionsCanonicalQ[exp, V[aa, {b, ib}, {c, ic}, d]];

swapInnerFermions@orderExtFields[orderF[orderB[a]]]

];

orderFermions[a___]:=(Message[orderFermions::syntax,a];Abort[])


(* Determine if inner fermion legs are in canonical order.
This means that the neighbours they connect to are ordered.
This canonical ordering is required to avoid wrong cancellations in the identification of diagrams that differ by a sign and the order of inner fermion legs. *)

(* do not apply to vertices with only two fermion legs *)
innerFermionsCanonicalQ[exp_op, v_V] /; 
  Count[v[[All, 1]], _?grassmannQ] < 4 := True
innerFermionsCanonicalQ[exp_op, v_V] := Module[
  {props, idRules, extFields, expId, extFields2ZeroRules, vConnectors,
    vNeighbours},
    
  (* get propagators and determine identification rules *)
  props = Cases[exp, P[___]];
  idRules = props /. P[a___] :> Rule @@ Sort[{a}];
  
  (* determine external fields (appear only once) *)
  extFields = 
   Select[Union@Cases[exp, {a_?fieldQ, b_}, Infinity], 
    Count[exp, #, \[Infinity]] == 1 &];
    
  (* rules to set all external fiels to zero *)
  
  extFields2ZeroRules = # :> Sequence[] & /@ extFields;
    
  (* expression with identifications, removed sf functions *)
  expId = exp /. P[___] :> Sequence[] /. idRules/. sf[_,_]:>Sequence[];
  
  (* get all internal fields of the vertex = connectors*)  
  vConnectors = List @@ (v /. extFields2ZeroRules /. idRules);
  
  (* get connected vertices *) 
  vNeighbours = 
   Function[conn, 
     Select[expId, Not@FreeQ[#, conn] && FreeQ[#, v] &]] /@ 
    vConnectors;
  
  (* delete internal fields and make lists of external legs of \
connected vertices *) 
  vNeighbours = 
   vNeighbours /. {a_?fieldQ, b_} :> 
       Sequence[] /; Not@MemberQ[extFields, {a, b}] /. op :> List /. 
    V :> Sequence;
  
  (* Is this list ordered? *)
  OrderedQ@vNeighbours
]


(* get the number of loops in a diagram *)

(* convert into a list, if sum is given *)
getLoopNumber[a_Plus]:=getLoopNumber/@List@@a;

(* get rid of numerical factors *)
getLoopNumber[a_Times]/;Count[a,op[__],\[Infinity]]==1:=getLoopNumber[Cases[a,op[__],\[Infinity]][[1]]];

getLoopNumber[a_op]:=Module[
	{props, vertices, regulatorInsertions},
	
	(* get number of vertices, propagators and regulator insertions *)
	vertices=Count[a, V[__] | S[__] | CO[__],\[Infinity]];
	props=Count[a, P[__],\[Infinity]];
	regulatorInsertions=Count[a, dR[__],\[Infinity]];

	(* number of loops *)
	props-vertices-regulatorInsertions+1
];

getLoopNumber[a___]:=Message[getLoopNumber::syntax,a];




(* ::Section:: *)
(* Derivation of a complete DSE *)


doDSE[a___,symmetry->broken,b___]/;FreeQ[{a,b},ansatz]:=(Message[doDSE::noAnsatz]; doDSE[a,symmetry->broken,b,ansatz->{}]);

(* if a real action is given transform it into a list *)
doDSE[action_,derivs_List,rest___,opts___?OptionQ]/;Cases[action,op[a__?(fieldQ@#[[0]](*check if arguments of op are fields*) &)],Infinity]=!={}:=
  doDSE[Replace[Union[Cases[{action},op[b___] :> Head /@ {b}, Infinity]],
  	{antiFermion_, fermion_} :> {fermion, antiFermion}(*order of fermions different in real and in symbolic action*),
  	 2],derivs,rest,opts];

(* if list of interactions given, create action first *)
doDSE[interactions_List,rest___,opts___?OptionQ]:=doDSE[generateAction[interactions],rest,opts];

(* only list of fields without indices, but $externalIndices is too short *)
doDSE[L_Times|L_Plus,derivs_List,b___]/;Depth@derivs==2&&Length@derivs>Length@$externalIndices:=Message[doDSE::tooFewExternalIndices,$externalIndices];

(* only list of fields without indices *)
doDSE[L_Times|L_Plus,derivs_List,b___]/;Depth@derivs==2:=doDSE[L,Transpose[{derivs,Take[$externalIndices,Length@derivs]}],b];

(* only list of list of fields without indices *)
doDSE[L_Times|L_Plus,derivs_List,b___]/;Depth@derivs==3&&And@@fieldQ/@Flatten@derivs:=doDSE[L,#,b]&/@derivs;

(* list of fields with indices *)
doDSE[L_Times|L_Plus,derivs_List,b___]/;Depth@derivs==4:=doDSE[L,#,b]&/@derivs;

doDSE[L_Times|L_Plus,derivs_List,b___]/;Or@@(fieldQ[#]&/@derivs[[All,2]]):=Message[doDSE::fieldAsIndex,Select[derivs[[All,2]],fieldQ[#]&]];

(* if no propagators are given, use only the ones given in the action *)
doDSE[L_Times|L_Plus,derivs_List,vertexTest_Symbol|vertexTest_Function,opts___?OptionQ]:=doDSE[L,derivs,
		Select[Cases[L, S[___], \[Infinity]], Length@# == 2 &]/. S[{Q1_, q1_}, {Q2_, q2_}] :> {Q1, Q2},vertexTest,opts];

doDSE[L_Times|L_Plus,derivs_List,allowedPropagators_List,vertexTest_Symbol|vertexTest_Function,opts___?OptionQ]:=Module[{firstDer,onePoint,multiPoint,
	sign,finalExp, complexFields},

(* get fields that are not necessarily fermions but directed, e.g., scalar complex fields;
this does not work if allowedPropagators is used (i.e. not {}) because then we cannot say what the
complex anti-field is  *)
complexFields={#,antiField@#}&/@Union@Cases[L, _?complexFieldQ,Infinity];

(* for 1PI vertex function add a minus sign due to its definition as the negative derivative of the effective action *)
sign= Which[Length@derivs>2,$signConvention (-1),True,(+1)];

(* first derivative *)
firstDer=deriv[L,First@derivs];

(* replace the fields and identify equal graphs *)
onePoint=replaceFields[firstDer];

(* perform additional differentiations and set sources to zero *)
multiPoint=deriv[onePoint,Sequence@@Rest@derivs];

(* get the correct sign due to fermions by ordering them; also order directed (complex) fields for identification*)
finalExp=If[sourcesZero/.Join[{opts},Options@doDSE],sortDummies@setSourcesZero[multiPoint,L,derivs,allowedPropagators,vertexTest,opts],
multiPoint,multiPoint];
finalExp=sign identifyGraphs[getSigns[finalExp], derivs];

finalExp

];

(* if no vertex test is given, use default *)
doDSE[L_Times|L_Plus,derivs_List,allowedPropagators___List,opts___?OptionQ]:=doDSE[L,derivs,allowedPropagators, (True&), opts]

doDSE[a___]:=Message[doDSE::syntax,a];




(* ::Section:: *)
(* Derivation of a complete RGE *)
(* approach using the logarithm *)


(* if a real action is given transform it into a list *)
doRGE[action_,derivs_List,rest___,opts___?OptionQ]/;Cases[action,op[a__?(fieldQ@#[[0]](*check if arguments of op are fields*) &)],Infinity]=!={}:=
  doRGE[Replace[Union[Cases[{action},op[b___] :> Head /@ {b}, Infinity]],
  	{antiFermion_, fermion_} :> {fermion, antiFermion}(*order of fermions different in real and in symbolic action*),
  	 2],derivs,rest,opts];

(* in the broken phase an explicit form for the average effective action has to be given *)
doRGE[interactions_List,derivs_List,rest___,opts___?OptionQ]:=Message[doRGE::noTruncation]/;Max[Length/@interactions]==2&&(symmetry/.{opts})===broken&&Not[Or@@(NumericQ/@interactions[[All,2]])];

(* if list of interactions given, create action first;
converting the list into an action is strictly speaking not required but 1) allows a more uniform approach and 2) allows to keep things parallel to DoDSE *)
doRGE[interactions_List,derivs_List,rest___,opts___?OptionQ]:=Module[{evenFields},
 evenFields=Cases[interactions,{Q_,even}:>Q];
 (* replace even fields definition by field two-point function *)
 doRGE[generateAction[interactions/.{Q_,even}:>{Q,Q}],derivs,rest,userEvenFields->evenFields,opts]
];

(* no derivatives, i.e., "zero-point function"; vertex test has no effect here *)
doRGE[L_Times|L_Plus,{},allowedPropagators_List,vertexTest_Symbol|vertexTest_Function,opts___?OptionQ]:=doRGE[L, {}, allowedPropagators, opts]

doRGE[L_Times|L_Plus,{},allowedPropagators_List,opts___?OptionQ]:=Module[{
	complexFields,zeroPoint,multiPointSources0,traceFieldRules,ind=insDummy[]},
	
(* get fields that are not necessarily fermions but directed, e.g., scalar complex fields;
this does not work if allowedPropagators is used (i.e. not {}) because then we cannot say what the
complex anti-field is  *)
complexFields={#,antiField@#}&/@Union@Cases[L, _?complexFieldQ,Infinity];

(* from the allowed propagators one can determine with which fields the trace can be done;
take twice the same field due to fermion ordering convention in propagators (since one index is in a propagator and one in a vertex) *)
traceFieldRules={{q1_,traceIndex1}:>{#,traceIndex1},{q2_,traceIndex2}:>{#,traceIndex2}}&/@Union@@allowedPropagators;

(* close the loop, ie also to bring in the minus sign from the supertrace *)
zeroPoint=Plus@@(1/2 op[dR[{$dummyField,traceIndex1},{$dummyField,ind}],P[{$dummyField,ind},{$dummyField,traceIndex2}]]/.traceFieldRules)
	/. P[Q1_, Q2_] :> -P[Q1, Q2] /; (Not@FreeQ[{Q1,Q2}, traceIndex1|traceIndex2,2] && (grassmannQ@Q1[[1]]||grassmannQ@Q2[[1]]))
	/.traceIndex1:>traceIndex2/.(Function[dField, 
   P[{dField[[2]], c_}, {dField[[1]], b_}] :> 
    P[{dField[[1]], b}, {dField[[2]], c}]] /@ Cases[complexFields,{_?(Head@#===boson&),_?(Head@#===boson&)}]);

(* order fermions and set sources to physical values *)
multiPointSources0= getSigns[setSourcesZeroRGE[zeroPoint,(*zeroSources,*)L,(*dirFields,*){{}},allowedPropagators,vertexTest,opts]];

identifyGraphs[sortDummies@multiPointSources0,{}]

];

(* only list of fields without indices, but $externalIndices is too short *)
doRGE[L_Times|L_Plus,derivs_List,b___]/;Depth@derivs==2&&Length@derivs>Length@$externalIndices:=Message[doRGE::tooFewExternalIndices,$externalIndices];

(* only list of fields without indices *)
doRGE[L_Times|L_Plus,derivs_List,b___]/;Depth@derivs==2&&derivs=!={}:=doRGE[L,Transpose[{derivs,Take[$externalIndices,Length@derivs]}],b];

(* only list of list of fields without indices *)
doRGE[L_Times|L_Plus,derivs_List,b___]/;Depth@derivs==3&&And@@fieldQ/@Flatten@derivs:=doRGE[L,#,b]&/@derivs;

(* list of fields with indices *)
doRGE[L_Times|L_Plus,derivs_List,b___]/;Depth@derivs==4:=doRGE[L,#,b]&/@derivs;


doRGE[L_Times|L_Plus,derivs_List,b___]/;Or@@(fieldQ[#]&/@derivs[[All,2]]):=Message[doRGE::fieldAsIndex,Select[derivs[[All,2]],fieldQ[#]&]];

(* if no propagators are given, use only the ones given in the action *)
doRGE[L_Times|L_Plus,derivs_List,vertexTest_Symbol|vertexTest_Function,opts___?OptionQ]:=doRGE[L,derivs,
		Select[Cases[L, S[___], \[Infinity]], Length@# == 2 &]/. S[{Q1_, q1_}, {Q2_, q2_}] :> {Q1, Q2},vertexTest,opts];


(* main code *)
doRGE[L_Times|L_Plus,derivs_List,allowedPropagators_List,vertexTest_Symbol|vertexTest_Function,opts___?OptionQ]:=Module[{
	onePoint,multiPoint,multiDer,sign,multiPointSources0,ind=insDummy[], extFields},

(* the external fields are given by the derivatives; extFields is just another name for it here *)
extFields=derivs;

(* for 1PI vertex function add a minus sign due to its definition as the negative derivative of the effective action *)
sign= Which[Length@derivs>2,$signConvention (-1),True,(+1)];

(* first derivative *)
onePoint=1/2 op[sf[First@derivs,{{$dummyField,traceIndex1},{$dummyField,ind}}],P[{$dummyField,traceIndex1},{$dummyField,ind}],-$signConvention V[First@derivs,{$dummyField,ind},{$dummyField,traceIndex2}]];

(* perform additional differentiations and set sources to zero *)
multiDer=derivRGE[onePoint,Sequence@@Rest@derivs];

(* if the sources are not set to 0, give back an intermediate results *)
multiPointSources0=sign If[sourcesZero/.Join[{opts},Options@doRGE],setSourcesZeroRGE[multiDer,L,extFields,allowedPropagators,vertexTest,opts],
multiDer,multiDer]/.traceIndex1:>traceIndex2;

(* closing the trace adds a minus sign if fields are fermionic; get other signs from fermions *)
multiPoint=getSigns[(multiPointSources0)/. P[Q1_, Q2_] :> -P[Q1, Q2] /; (Not@FreeQ[{Q1,Q2}, traceIndex1|traceIndex2,2] && (grassmannQ@Q1[[1]]||grassmannQ@Q2[[1]]))];

multiPoint=sortCanonical[multiPoint, derivs];

(* the first dummy sorting is required for the identification to work; the second one treats the dummies introduced by derivPropagatorsdt *)
If[tDerivative/.Join[{opts},Options@doRGE],
	sortDummies[identifyGraphs[multiPoint,extFields]/.a_op:>derivPropagatorsdt@a],
	sortDummies[identifyGraphs[multiPoint,extFields]],
	sortDummies[identifyGraphs[multiPoint,extFields]/.a_op:>derivPropagatorsdt@a]
]

];

(* no vertex test given, use default *)
doRGE[L_Times|L_Plus,derivs_List,allowedPropagators___List,opts___?OptionQ]:=doRGE[L, derivs, allowedPropagators, (True&), opts]

doRGE[a___]:=Message[doRGE::syntax,a];


(* derivation of correlation functions for composite operators *)

(* main code *)
doCO[action_Times|action_Plus|action_List, compOp_, filter_:(True&), opts___?OptionQ]:=Module[{
	singleExp, allFields, extFields, compOpFieldsRep, compOpFieldsRepS0, compOpFieldsRepS0Filtered(*, compOpFieldsRepS0FilteredId*)
	},
	
	singleExp = Expand[compOp];
	singleExp = If[ListQ@singleExp||Head[singleExp]==Plus, singleExp[[1]], compOp];
	allFields = Cases[singleExp, {_?fieldQ, _}, Infinity];
	extFields = Select[allFields, Count[allFields, #]==1&];
		
	compOpFieldsRep = replaceFields@compOp;
	compOpFieldsRepS0 = setSourcesZero[compOpFieldsRep, action, {}]	;
	compOpFieldsRepS0Filtered = Select[compOpFieldsRepS0, filter];
	(* currently the identification does not work reliably for disconnected diagrams, thus deactivate *) 
	(* compOpFieldsRepS0FilteredId = identifyGraphs[compOpFieldsRepS0Filtered, extFields];*)
	compOpFieldsRepS0Filtered	
]

doCO[a___]:=Message[doCO::syntax,a];




(* ::Section:: *)
(* Control Functions *)


(* check if no indices appear more often than twice;
indicesTest is the test function and checkIndices performs the test *)

indicesTest[a_]/;Not@FreeQ[a,{_?fieldQ,_}]:=Module[
{inds},

inds=Transpose[Cases[a,{_,ind_},Infinity]][[2]];
Union@Select[{#,Count[inds,#]}&/@inds,#[[2]]>2&][[All,1]]
];

indicesTest[a_]:={};


checkIndices[a_op]:=checkIndices[{a}];

checkIndices[a_]/;(b=Select[Cases[a,op[___],\[Infinity]],(indicesTest[#]!={}&)])!={}:=
	Message[checkIndices::multipleIndices,b,indicesTest[b]];

checkIndices[a_]:=Message[checkIndices::ok];



(* check if the syntax is correct, i.e. all elements of op are propagators, vertices or fields,
the propagators have the correct syntax of two indices and the vertices that of several indices *)


opTest[a_]:=Select[Cases[{a},op[___],\[Infinity]],Cases[#,_S|_V|_P|_List]!= List@@#&];


propagatorTest[a_]:=Select[Cases[{a},P[___],\[Infinity]],Not@MatchQ[#,P[{_?fieldQ,_},{_?fieldQ,_}]]&];


vertexTest[a_]:=Select[Cases[{a},S[___]|V[___],\[Infinity]],Cases[#,_List]!= List@@#&];


regulatorInsertionTest[a_]:=Select[Cases[{a},dR[___],\[Infinity]],Not@MatchQ[#,dR[{_?fieldQ,_},{_?fieldQ,_}]]&];


checkSyntax[a_op]:=checkSyntax[{a}];

checkSyntax[a_]/;(b=opTest[a])!={}:=Message[checkSyntax::op,b];

checkSyntax[a_]/;(b=propagatorTest[a])!={}:=Message[checkSyntax::propagator,b];

checkSyntax[a_]/;(b=vertexTest[a])!={}:=Message[checkSyntax::vertex,b];

checkSyntax[a_]/;(b=regulatorInsertionTest[a])!={}:=Message[checkSyntax::regulatorInsertion,b];

checkSyntax[a_]:=Message[checkSyntax::ok];



(* all checks in one function *)

checkAll[a_]:=(checkIndices[a];checkSyntax[a];checkFields[a];);


(* checks in the action:
all indices have to appear twice (no free indices) *)

indicesTwice[a_]:=Module[{inds},
inds=Cases[a,{b_,c_}:>c,\[Infinity]];
Complement[inds,Select[inds,Count[inds,#]==2&]]
];


checkAction[a_op]:=checkAction[{a}];

checkAction[a_]/;(b=Select[Cases[a,op[___],\[Infinity]],indicesTwice[#]!={}&])!={}:=Message[checkAction::error,b,indicesTwice[b]];

checkAction[a_]:=(checkSyntax[a];checkFields[a];Message[checkAction::ok]);


checkFields[a_op]:=checkFields[{a}];

checkFields[a_]/;(d=Union@Select[Cases[a,{b_,c_}:>b,\[Infinity]],Not@fieldQ[#]&])!={}:=Message[checkFields::undefinedField,d];

checkFields[a_]:=Message[checkFields::ok];




(* ::Section:: *)
(* Functions for Output *)


(* writing the output with shorter indices so that it looks nicer *)


shortExpressionSingle[S[a__]]:=Subsuperscript[$bareVertexSymbol,StringJoin[ToString/@Riffle[{a}[[All,1]]," "]],StringJoin[ToString/@Riffle[{a}[[All,2]]," "]]];

shortExpressionSingle[{Q_,q_}]:=Subscript[Q,q];

shortExpressionSingle[V[a__]]:=Subsuperscript[$vertexSymbol,StringJoin[ToString/@Riffle[{a}[[All,1]]," "]],StringJoin[ToString/@Riffle[{a}[[All,2]]," "]]];

shortExpressionSingle[CO[a__]]:=Subsuperscript[$compOpSymbol,StringJoin[ToString/@Riffle[{a}[[All,1]]," "]],StringJoin[ToString/@Riffle[{a}[[All,2]]," "]]];

shortExpressionSingle[dR[a__]]:=DisplayForm@RowBox[{SubscriptBox["\[PartialD]", "t"],Subsuperscript[$regulatorInsertionSymbol,StringJoin[ToString/@Riffle[{a}[[All,1]]," "]],StringJoin[ToString/@Riffle[{a}[[All,2]]," "]]]}];

shortExpressionSingle[P[{F1_,i1_},{F2_,i2_}]]:=Subsuperscript[$propagatorSymbol,ToString@F1<>" "<>ToString@F2,ToString[i1]<>" "<>ToString[i2]];

shortExpressionSingle[sf[_,_]] := 1


shortExpressionStyle[a_,opts___]:=Style[a,
		Join[{opts}/.{Rule[_,_]:>Sequence[],RuleDelayed[_,_]:>Sequence[]},
			FilterRules[Join[Cases[{opts},Rule[_,_]|RuleDelayed[_,_]],Options[shortExpression]],Options@Style]]];


(*shortExpression[a_Times,opts___]:=shortExpressionStyle[shortExpression[#,opts]&/@a,opts];*)
(* handles signs correctly *)
shortExpression[a_Times,opts___]:=a[[1]] shortExpression[a[[2]],opts];

shortExpression[a_Plus,opts___]:=shortExpressionStyle[shortExpression[#,opts]&/@a,opts];

shortExpression[a_?NumericQ,opts___]:=shortExpressionStyle[a,opts];

shortExpression[op[a__],opts___]:=shortExpressionStyle[Times@@(shortExpressionSingle/@(op[a]/.op:>List)),opts];

shortExpression[a___]:=Message[shortExpression::syntax,a];


(* identification of diagrams via their names type *)
(* Extract number and multiplicities of vertices. *)
getVertexNumbers[n_?NumericQ a_] := getVertexNumbers[a]

getVertexNumbers[a_op] := Sort@Cases[a, V[v__]|S[v__] :> Length@{v}, Infinity]


(* Functions to determine the diagram type. *)

getDiagramType[n_?NumericQ a_] := getDiagramType[a]

getDiagramType[a_op] := diagramTypes[{getLoopNumber[a], getVertexNumbers[a]}]


(* Extract a certain diagram type. *)

extractDiagramType[diags_, type_String] := Select[diags, getDiagramType[#] == type &]


(* Group diagrams by classes. *)

groupDiagrams[exp_] := (Function[class, {class, extractDiagramType[exp, class]}] /@ Values@diagramTypes) /. {_, 0} :> Sequence[]


(* vertexDummies: auxiliary functions for DSEPlotList and RGEPlotList, determines unique dummies for the vertices
so that they can be used in DSEPlotList/RGEPlotList as points *)

vertexDummies[a_?NumericQ b_,opts___?OptionQ] := {vertexDummies[b,opts], a};

vertexDummies[a_Plus,opts___?OptionQ] := vertexDummies[#,opts] & /@ List @@ a;

(* the function sort allows to give a function applied on the list of vertices;
up to now only necessary when invoked from DSEPlotCompare because a unique vertex representative is needed;
otherwise use an "empty" function *)
(*vertexDummies[a_op,fields_List,sort_:(#&),opts___?OptionQ]*)
vertexDummies[a_op,sort_:(#&),opts___?OptionQ] /;Not@OptionQ@sort:=(*vertexDummies[a,opts]=*) Module[
  {allIndices, indices, vertices, bareVertexRepres, CORepres, pureProps, purePropRepres, bareVertexRule, CORule, purePropRule,
  	externalFields, externalFieldsSingle, externalPropagators, externalPropagatorsSingleFields, propagators, legs, regulators,regulatorsRepres,regulatorRule,
  	verticesRepres, verticesRepresList, identificationList,sameVertTest, extText, extFieldText, allDirFields},
  
  allDirFields=getDirectedFieldsList[a];

  (* don't plot the index when called from DSEPlotCompare because then graphs won't be identified *)
  extText[q_]:=" ext "<>ToString@q;
  extText[q_]/;sort===Sort:=" ext ";
  extFieldText[q_]:=" exField "<>ToString@q;
  extFieldText[q_]/;sort===Sort:=" extField ";

  (* get all vertices from the epxression and standalone propagators *)
  vertices = Cases[a, _V | _S | _dR | _CO | _List];

  bareVertexRepres=Cases[a,_S][[All,1]];
  
  CORepres = Cases[a, _CO][[All,1]];
  
  (* get all indices and identify those at the same vertex *)
  allIndices = Cases[a, {_?fieldQ, _}, Infinity];
  sameVertTest[c_,d_]:=(Or @@ Function[{vert}, Not@FreeQ[vert, c] && Not@FreeQ[vert, d]] /@ vertices );
  indices = Union[allIndices, SameTest -> sameVertTest ];
  
  (* get legs *)
  legs = Select[allIndices, Count[a, #, Infinity] == 1 &];

  (* regulator insertions *)
  regulators=Cases[a, _dR ];
  regulatorsRepres=Cases[a,_dR][[All,1]];
  regulatorRule=#->(#/.{Q_,q_}:> {Q,q,"dt R"})&/@regulatorsRepres;

  (* external fields *)
  externalFields = Cases[a, {_?fieldQ, _}];
  
  (* external fields not coupled to anything *)
  externalFieldsSingle = Intersection[externalFields, legs];
  externalFields = Complement[externalFields, externalFieldsSingle];
  
  (* delete exernal fields from legs *) 
  legs = Complement[legs, externalFieldsSingle];
  
  (* get all internal propagators; sort according to order of directed fields *)
  propagators =Sort[#,MemberQ[allDirFields,{#1,#2}]&]&/@ Cases[a, _P];
    
  (* pure propagators *)
  pureProps = Cases[propagators, P[{f1_,i1_},{f2_,i2_}]/;Not@FreeQ[legs,i1]&&Not@FreeQ[legs,i2]];
  
  (* add pure propagators to vertices *)
  vertices = Join[vertices, pureProps];
  purePropRepres = pureProps[[All,1]];
  
  (* drop standalone propagators *)
  propagators = Complement[propagators, pureProps];

  (* add the propagators to external points *)
  externalPropagators = Join[legs /. {Q_, q_?(Not@ListQ@# &)} :> P[{Q, q}, {antiField@Q, q, " leg "<>ToString@q}], externalFields /. {Q_, q_?(Not@ListQ@# &)} :> P[{Q, q}, {antiField@Q, q,extText[q]}]];
  
  (* for single external fields create a propagator which will be plotted transparent *)
  externalPropagatorsSingleFields = externalFieldsSingle /. {Q_, q_?(Not@ListQ@# &)} :> P[{Q, q}, {antiField@Q, q,extFieldText[q]}];

  propagators =Sort[#,(* this is indeed the correct direction for directed fields *)Not@MemberQ[allDirFields,{#1[[1]],#2[[1]]}]&]&/@ Join[propagators, externalPropagators, externalPropagatorsSingleFields];

  (* give all propagators a "name" *)
  propagators=propagators/.P[{B_, b_,bl___}, {C_, c_,cl___}]:>P[{B, b,bl}, {C, c,cl},StringJoin[ToString@B,bl," ",ToString@C,cl]];

  (* get indices of all vertices and choose one representative;
     sort to have uniquely defined representatives when plotting with DSEPlotCompare, but not otherwise since then the bare vertex can "disappear" *)
  verticesRepres = {First@#, List @@ #} & /@ sort/@ vertices;
  (* determine how to add an "S", "CO" or "P" to the bare vertex representative; this is used in the VertexRenderingFunction to plot S, P and V differently *)
  bareVertexRule=#->(#/.{Q_,q_}:> {Q,q,"S"})&/@bareVertexRepres;
  CORule=#->(#/.{Q_,q_}:> {Q,q,"CO"})&/@CORepres;
  purePropRule = #->(#/.{Q_,q_}:> {Q,q,"P"})&/@purePropRepres;

  (* prepare a list for the identification of every index with the representative; Hold necessary for the use of Riffle *)
  verticesRepresList = Flatten[Partition[Riffle[verticesRepres[[#, 2]], Hold@verticesRepres[[#, 1]], {2, -1, 2}], 2] & /@ Range@Length@verticesRepres // ReleaseHold, 1];

  (* identification rules *)
  identificationList = Thread[Rule[#[[1]], #[[2]]], 1] & /@ verticesRepresList;

  (* use propagators for the rules of the Graph *)
  propagators=propagators /. identificationList /.bareVertexRule/.CORule/.purePropRule/.regulatorRule/. P[{B_, b_, bl_: ""}, {C_, c_, cl_: ""},d_] :> {Rule[{B, b, bl}, {C, c, cl}],d}
   
  ];
  

(* auxiliary function for DSEPlotList for plotting arrows/lines for fermions/bosons instead of lines *)

arrowLine[field_, coords___, opts___?OptionQ]:=Module[{},
	
	If[grassmannQ[field] || complexFieldQ[field] || antiComplexFieldQ[field],
		GraphElementData["Arrow"][coords], 
    	GraphElementData["Line"][coords]
    ]
];



(* auxiliary function for comparison of graphs; derived from DSEPlotList *)
Clear@DSEPlotCompare;
DSEPlotCompare[a_op]:=DSEPlotCompare[a]=GraphPlot[vertexDummies[a,Sort]];


(* overview over DSEPlot functions:
DSEPlotCompare: for comparing Graphs (deprecated)
DSEPlot: main command, plots complete DSEs
DSEPlotList: "old" main command, does the plotting, but gives a list
DSEPlotGrid: puts the graphs into a GridBox
*)


(* auxiliary function for Graph: sets the styles of edges for fields *)

(* determines the style of an edge based on given style rules and \
grassmann/non-grassmann *)
getEdgeShapeFunction[edge_, plotStyles_List] := Module[{field},
  
  (* extract which field we have in the propagator; 
  make detour via the label, because the vertices contain one representative field only;
  note that for mixed propagators only one field is taken (could be taken into account, though, by writing a dedicated plot function *)
  field = ToExpression@
    StringCases[{edge[[2]]}, 
      f__ ~~ " " ~~ __ /; fieldQ[ToExpression[f]] :> f][[1, 1]];

  (* extract the proper style and set the arrow *)
  Join[Cases[plotStyles, {field | antiField@field, 
      style___} :> Unevaluated[({style, arrowLine[field, ##]} &)]], {({arrowLine[field, ##]}&)}(* in case no plotStyle matches *)][[1]]
]


(* auxiliary function for Graph: sets the styles of vertices *)

getVertexShapeFunction[{x_,y_}, label_, {w_, h_}, opts___?OptionQ]/;Not@StringFreeQ[label[[3]],"S"] := 
	bareVertexSymbol[{x,y},{w,h}]/.Join[{opts},Options@DSEPlotList]

getVertexShapeFunction[{x_,y_}, label_, {w_, h_}, opts___?OptionQ]/;Not@StringFreeQ[label[[3]],"dt R"] :=
	regulatorSymbol[{x,y},{w,h}]/.Join[{opts},Options@DSEPlotList]

getVertexShapeFunction[{x_,y_}, label_, {w_, h_}, opts___?OptionQ]/;Not@StringFreeQ[label[[3]],"CO"] :=
	coSymbol[{x,y},{w,h}]/.Join[{opts},Options@DSEPlotList]

getVertexShapeFunction[{x_,y_}, label_, {w_, h_}, opts___?OptionQ]/;Not@StringFreeQ[label[[3]],"leg"] :=
	{Text[Style[StringReplace[label[[3]],"leg":> ""],Sequence@@(indexStyle/.Join[{opts},Options@DSEPlotList]),FontSize:>14],{x,y}+{0,0.3}]}

getVertexShapeFunction[{x_,y_}, label_, {w_, h_}, opts___?OptionQ]/;Not@StringFreeQ[label[[3]],"ext"] :=
	{Text[Style[StringReplace[label[[3]],"ext":> ""],Sequence@@(indexStyle/.Join[{opts},Options@DSEPlotList]),FontSize:>14],{x,y}+{0,0.3}],Circle[{x,y},0.1]}

getVertexShapeFunction[{x_,y_}, label_, {w_, h_}, opts___?OptionQ] := vertexSymbol[{x,y},{w,h}]/.Join[{opts},Options@DSEPlotList]


(* Plotting graphs using Graph; if a list of directed propagators/fermions is given, arrows will be used;
the plotting is done in DSEPlotList, which creates a list of plots
the user invokes DSEPlot and gets a complete DSE, but with the option output -> List one can also get a list *)


(* with edges rendered specially *)
(* if action is given as list, it has to be transformed first *)
DSEPlotList[a_List,plotRules_List,opts___?OptionQ]/;And @@ (fieldQ /@ Flatten[a]):=
	DSEPlotList[generateAction[a],plotRules,opts];

DSEPlotList[a_,plotRules_List,opts___?OptionQ]/;FreeQ[a,Rule,Infinity]:=
	DSEPlotList[vertexDummies[a,opts],plotRules,opts];

DSEPlotList[{a_List,b_?NumericQ},plotRules_List,opts___?OptionQ]:=Module[{graph, exponent, verts,
	plotRulesAll, vsize, edgeLabels, multiEdges, multiEdgeRules, multiEdgeReps},

	(* sort real complex fields for easier identification below, take the field info from the label *)
	graph = Replace[a, {c_Rule, d_}:>{Sort@c, d}/;bosonQ[ToExpression@StringCases[{d}, f__ ~~ " " ~~ __ /; fieldQ[ToExpression[f]] :> f][[1, 1]]], {1}];
	
	(* In Mathematica 12.0.0 EdgeRenderingFunction was superseded by EdgeShapeFunction at it was modified getting less arguments.
	The label is no longer provided as argument and thus edges can no longer be distinguished properly.
	To avoid problems in future versions, EdgeShapeFunction is used together with Graph (instead of GraphPlot).
	The remaining problem is that edges connecting the same vertice cannot be distinguished. Possibly this is a bug/unforeseen problem in M12.
	In case of two different edges a cheat is to use UndirectedEdge and DirectedEdge, which can be distinguished (the latter even twice because of two directions;
	however, this is not used as different directions can appear by themselves).
	If there are more different edges, no method is known. A warning is given then. *)
	
	(* get all adges that connect the same vertices *)
	multiEdges = Select[graph, Count[graph, #[[1]], \[Infinity]] >= 2 &];
	(* group them *)
	multiEdges = Union/@GatherBy[multiEdges, #[[1]]&];
	
	(* make sure only one representative is there *)
	multiEdgeReps = Function[d, DeleteDuplicatesBy[d, #[[2]]&]]/@multiEdges;
	
	(* give a warning if there are more than 2 different edges *)
	If[Or@@(Length[#]>2 & /@ multiEdgeReps), Message[DSEPlotList::multiPropagators]];
	
	(* rules for replacing some of the edges: use undirected edges *)
	multiEdgeRules = (#->{UndirectedEdge@@#[[1]], #[[2]]})&/@Select[multiEdgeReps, Length@#>1&][[All,2]];

	(* replace the edges *)
	graph = graph/.multiEdgeRules;
	
	(* add field style to propagators; if no style is given, add edge labels *)
	If[plotRules==={},
		edgeLabels = (#[[1]]->#[[2]])&/@graph,
		edgeLabels = {}
	];
	
	(* set properties of edges *)
	graph = Property[#[[1]], EdgeShapeFunction -> getEdgeShapeFunction[#, plotRules]] & /@ graph;
	
	(* add vertex style *)
	verts = Union@@List@@@a[[All,1]];
	graph = {Function[vert, Property[vert, VertexShapeFunction->(getVertexShapeFunction[##, opts]&)]] /@ verts, graph};
	
	(* extend plot Rules also to antifields *)
	plotRulesAll=Union@Replace[plotRules, {c_?fieldQ, d__} :> Sequence[{c, d}, {antiField@c, d}], 1];
	
	(* get fields that are not necessarily fermions but directed, e.g., scalar complex fields *)
	(*dirFields=(directedFields/.Join[{opts},Options@DSEPlot]);*)
	
	(* a smaller regulator symbol size is required for zero leg graphs in RGEs to avoid that the regulator symbol is larger than the loop *)
	vsize = 0.15;
	If[Not@FreeQ[a,"dt R"] && And@@(StringFreeQ[#, "leg"]&/@a[[All,2]]), vsize=0.03];
	
	(* exponent -1 for inverse propagators *)
	exponent=If[Length@a===2 && FreeQ[a, "P"],"-1","",""];
	
	Labeled[
	Graph[Sequence@@graph, 
			FilterRules[Join[{opts},Options@DSEPlot],Options@Graph],
			VertexSize->vsize,
			GraphLayout -> "SpringElectricalEmbedding",
			EdgeLabels -> edgeLabels],
			(* for positive integers explicitly print the +, for positive Rationals also, but it has to be prevented that the + goes into the numerator;
				RowBox necessary to prevent automatic ordering *)
	{Style[b/.{
	 1 :> Style["+", ShowStringCharacters -> False],
	 (c_Integer | c_Real | c_Rational)?Positive :> DisplayForm@RowBox[{Style["+", ShowStringCharacters -> False], c}],
	 -1 :> Style["-", ShowStringCharacters -> False]},(factorStyle/.Join[{opts},Options@DSEPlot])],
	 Overscript[Style["",FontSize:>50,ShowStringCharacters -> False],Style[exponent,ShowStringCharacters -> False,(factorStyle/.Join[{opts},Options@DSEPlot])]]},
	 {Left,Right}
	 ]
];

DSEPlotList[a_,plotRules_List,opts___?OptionQ]/;FreeQ[a,Rule[_,_],2]:=DSEPlotList[#,plotRules,opts]&/@a;

DSEPlotList[a_List,plotRules_List,opts___?OptionQ]:=DSEPlotList[{a,1},plotRules,opts];

(* without edge rendering *)
DSEPlotList[a_, opts___?OptionQ]:=DSEPlotList[a,{},opts];

DSEPlotList[a___]:=Message[DSEPlot::syntax,a];



(* plots are all based on DSEPlot, they are only distinguished by the option type *)

COPlot[args___]:=DSEPlot[args, type->"CO"];

RGEPlot[args___]:=DSEPlot[args, type->"RGE"];


(* plot the complete equation including the left-hand side; employ a grid;
if no PlotRules are given, call DSEPlotList accordingly without it *)

(* if there is only a number *)
DSEPlot[a_?NumericQ,___]:=a;

(* fields not defined *)
DSEPlot[a_,rest___]/;Not[And@@(fieldQ/@Union[Cases[a,{b_,_}:>b,{2,Infinity}]])]:=(Message[DSEPlot::fieldsUndefined,a];Abort[]);

(* if only style definitions for one field are given and one pair of brackets is missing rewrite it with the proper syntax *) 
DSEPlot[a_,{f_?fieldQ,styleDefs___},rest___]:=DSEPlot[a,{{f,styleDefs}},rest];

(* if plotRules are not properly defined *)
DSEPlot[a_,plotRules_List,___]/;Not@MatchQ[plotRules, {{__},___}]:=Message[DSEPlot::plotRules,plotRules];

DSEPlot[a_,plotRules_List:{},len_Integer:5,opts___?OptionQ]/;(output/.Join[{opts},Options@DSEPlot])===List:=
	DSEPlotList[a,plotRules/.{}:>Sequence[],opts];

(* plot lists of diagrams as such *)
DSEPlot[a_List,plotRules_List:{},len_Integer:5,opts___?OptionQ]:=DSEPlot[#,plotRules,len,opts]&/@a;

(* plot sum of op operators; normally a single graph is plotted alone, except the option output is set to forceEquation *)
DSEPlot[a_,plotRules_List:{},len_Integer:5,opts___?OptionQ]/;And@@(Not@FreeQ[#, op[__]] & /@ List@@Expand[a])||(output/.Join[{opts},Options@DSEPlot])===forceEquation:=Module[
	{expandeda, rhs, lhs, inds, lhsFields, exponent, Q,q, eqType, lhsSymbol},
	
	(* what equation should be plotted? *)
	eqType = type/.Join[{opts}, Options[DSEPlot]];
	
	(* expand a or otherwise there may be problems with parentheses *)
	expandeda=Expand@a;
	
	(* get all fields and indices *)
	inds=Cases[First@Cases[{expandeda}, op[___], \[Infinity]], {Q_, q_}, \[Infinity]];
	
	(* determine the externals *)
	lhsFields=Select[inds, Count[inds, #] == 1 &];
	
	(* symbol to use for lhs *)
	lhsSymbol = If[eqType==="CO",
		CO,
		V,
		V];
	
	(* plot the left-hand side; take only the graph without prefactor; if the diagram has no external legs, take just Gamma_k; not required for DSEs as there are no vacuum diagrams *)
	lhs=Which[lhsFields=={},
		DisplayForm@StyleBox[SuperscriptBox["\[CapitalGamma]", "k"],factorStyle/.Join[{opts},Options@RGEPlot],FontSize->20],
		True,
		DSEPlotList[op@lhsSymbol[Sequence@@lhsFields], plotRules/.{}:>Sequence[],opts][[1]]
	];
	
    (* exponent -1 for propagators *)
    exponent=If[Length@lhsFields===2,"-1","",""];

	(* plot the right hand side *)
	rhs=DSEPlotList[expandeda,plotRules/.{}:>Sequence[],opts];

	Which[eqType == "DSE",
		DSEPlotGrid[rhs,{lhs,exponent},len,opts],
		eqType == "RGE",
		RGEPlotGrid[rhs,{lhs,exponent},len,opts],
		eqType == "CO",
		DSEPlotGrid[rhs,{lhs,exponent},len,opts]
		]
];

(* plot single diagrams as such *)	
DSEPlot[a_,plotRules_List:{},len_Integer:5,opts___?OptionQ]/;Count[a, op[___], \[Infinity]]==1||Head@a==op:=DSEPlotList[a,	plotRules/.{}:>Sequence[],opts];

DSEPlot[a___]:=Message[DSEPlot::syntax,a];


(* for DSEs and CO *)
DSEPlotGrid[rhs_,  {lhs_,exponent_}, len_, opts___?OptionQ] := 
  Module[{partitioned, lhsLabeled, i},

   (* divide into the correct length; if rhs is only one graph convert it to a list *)  
    partitioned =  Insert[Partition[Flatten[{rhs}], len - 1, len - 1, 1, Style["",ShowStringCharacters->False]], Style["",ShowStringCharacters->False], 
     Table[{i, 1}, {i, Ceiling[Length@rhs/(len - 1)]}]];

   (* create the equal sign as label for lhs and add -1 exponent for propagator *)
   lhsLabeled = Labeled[lhs, Style[Row[{Overscript[Style["",FontSize:>50],Style[exponent,(factorStyle/.Join[{opts},Options@DSEPlot])
   			/.(FontSize:>w_Integer):>(FontSize:>w)]],"="}],
   		ShowStringCharacters->False,factorStyle/.Join[{opts},Options@DSEPlot]], Right];

   (* put in the lhs *)
   ReplacePart[Grid[Sequence @@@ List /@ partitioned,FilterRules[{opts},Options@Grid]], 
    lhsLabeled, {1, 1, 1}]
   
];


RGEPlotGrid[rhs_,  {lhs_,exponent_}, len_, opts___?OptionQ] := 
  Module[{partitioned, lhsLabeled, i},

   (* divide into the correct length; if rhs is only one graph convert it to a list *)  
    partitioned =  Insert[Partition[Flatten[{rhs}], len - 1, len - 1, 1, Style["",ShowStringCharacters->False]], Style["",ShowStringCharacters->False], 
     Table[{i, 1}, {i, Ceiling[Length@rhs/(len - 1)]}]];

   (* create the equal sign as label for lhs and add -1 exponent for propagator *)
   lhsLabeled = Labeled[lhs,
   			Style[#,ShowStringCharacters->False,factorStyle/.Join[{opts},Options@DSEPlot]]&/@
   				{DisplayForm@SubscriptBox["\[PartialD]", "t"],Row[{Overscript[Style["",FontSize:>50],Style[exponent,(factorStyle/.Join[{opts},Options@DSEPlot])
   			/.(FontSize:>w_Integer):>(FontSize:>w)]],"="}]}, {Left,Right}];

   (* put in the lhs *)
   ReplacePart[Grid[Sequence @@@ List /@ partitioned,FilterRules[{opts},Options@Grid]], 
    lhsLabeled, {1, 1, 1}]
   
];


(* symbols for plotting *)

boxSymbol[x_,{w_,h_}]:={GrayLevel[0.4],Rectangle[x-{w, h},x+{w, h}]};
boxSymbol[x_]:=boxSymbol[x, {1,1}];

crossSymbol[x_,{w_,h_}]:=Module[{rad=Sqrt[w^2+h^2]},
	{
	Circle[x, rad], 
    Line[{x+rad{-Sin[\[Pi]/4], -Cos[\[Pi]/4]}, x+rad{Sin[\[Pi]/4], Cos[\[Pi]/4]}}], 
    Line[{x+rad{Sin[\[Pi]/4], -Cos[\[Pi]/4]},x+rad {-Sin[\[Pi]/4], Cos[\[Pi]/4]}}]
	}
];
crossSymbol[x_]:=crossSymbol[x,{1,1}];

diskSymbol[x_,{w_,h_}]:={Disk[x,Sqrt[w^2+h^2]]};
diskSymbol[x_]:=diskSymbol[x,{1,1}];

diskOpenSymbol[x_,{w_,h_}]:={Circle[x, Sqrt[w^2+h^2]]};
diskOpenSymbol[x_]:=diskOpenSymbol[x, {1,1}];

diskTinySymbol[x_,{w_,h_}]:={Disk[x,0.2*Sqrt[w^2+h^2]]};
diskTinySymbol[x_]:=diskTinySymbol[x,{1,1}];

triangleSymbol[x_,{w_,h_}]:={GrayLevel[0.4],Polygon[{x- {w, h}, x + {w, 0}, x + {-w, h}}]};
triangleSymbol[x_]:=triangleSymbol[x,{1,1}];


(* choice of different regulators *)

(* these names are kep for compatibility *)

regulatorBox[x___]:=boxSymbol[x];
regulatorCross[x___]:=crossSymbol[x];




(* ::Section:: *)
(* Tools *)


(* create the propagator rules for the dummy field according to the allowed propagators given in propagatorList
in the form {{A,A},{c,cb}} *)
createPropagatorRules[propagatorList_List, fields_List,opts___?OptionQ]:=Module[{ffields,propagators,propsClasses},

ffields=Flatten@fields;

(* create all possible propagators, i.e. when {c,cb} is given also create {cb,c}*)
propagators=P[{#[[1]],ind1$}, {#[[2]],ind2$}] & /@ Union[propagatorList,Reverse/@propagatorList];

(* classes: a field and all the propagators where it is involved *)
propsClasses={#,Cases[propagators, P[{#,_},{_?fieldQ,_}]]}&/@ffields;

(* all possible propagator replacements that are allowed *)

Join[ 
	Function[{Q},P[{Q[[1]],ind1_},{$dummyField,ind2_}]->Evaluate@Apply[Plus,Q[[2]]]]/@propsClasses,
	Function[{Q},P[{$dummyField,ind2_},{Q[[1]],ind1_}]->Evaluate@Apply[Plus,Reverse/@Q[[2]]]]/@propsClasses,
	{P[{$dummyField,ind1_},{$dummyField,ind2_}]->Evaluate[Plus@@(propagators/.P[{Q1_,Q2_}]:>P[{Q1,ind1},{Q2,ind2}])]}
]

];



(* get all fields which are directed, i.e., fermions and complex bosonic fields *)
getDirectedFieldsList[exp_]:=Module[{fermions,antiFermions,complexFields},
	fermions=Cases[exp,{a_?fermionQ,___}:>{a,antiField[a]},Infinity];
	(* need to take into account also anti-fermions as the order in V and P differs *)
	antiFermions=Cases[exp,{a_?antiFermionQ,___}:>{antiField[a],a},Infinity];
	complexFields=Cases[exp,{a_?complexFieldQ,___}:>{a,antiField[a]},Infinity];
	Union[fermions,antiFermions,complexFields]
];


(* Fields: Properties of fields need to be defined.
The following types of fields are defined:
Real bosons, fermions (or Grassmann fields in general), complex bosons.
Fields must be defined before doing any calculations.
No other function should interfer with the fields' definitions. *)

(* error if the dummy field would be overwritten *)
setFields[fields___] /; Not[FreeQ[{fields}, $dummyField]] := Message[setFields::dummyField]

(* some defaults *)
setFields[bosons_List] := setFields[bosons, {}, {}]
setFields[bosons_List, fermions_List] := 
 setFields[bosons, fermions, {}]
(* check syntax of fermionic and complex fields (need to be pairs) *)

setFields[bosons_List, fermions_List, complexFields_List] /; 
   Not[And @@ 
     Flatten[MatchQ[#, {_, _}] & /@ Join[fermions, complexFields]]] :=
   Message[setFields::noPair, fermions, complexFields];
setFields[bosons_List, fermions_List, complexFields_List] := Module[{},
  
  (* set the field types *)
  (* difference to DoFun2: 
  complex fields are their own types *)
  (# /: fieldType[#] := boson) & /@ bosons;
  (# /: fieldType[#] := fermion) & /@ fermions[[All, 1]];
  (# /: fieldType[#] := antiFermion) & /@ fermions[[All, 2]];
  (# /: fieldType[#] := complex) & /@ complexFields[[All, 1]];
  (# /: fieldType[#] := antiComplex) & /@ complexFields[[All, 2]];
 
  (* set all fields to Head field: Not done, 
  since this does not fully work as expected *)
  (*(#/:Head[#]=
  field)&/@Flatten[{bosons,fermions,complexFields}];*)
  
  (* set (anti-)commutating property *)
  (# /: grassmannQ[#] = False) & /@ Flatten[{bosons, complexFields}];
  (# /: grassmannQ[#] = True) & /@ Flatten[fermions];
  (# /: cFieldQ[#] = True) & /@ Flatten[{bosons, complexFields}];
  (# /: cFieldQ[#] = False) & /@ Flatten[fermions];
  
  (* define the anti-  fields *)
  (antiField[#[[1]]] = #[[2]]) & /@ Transpose[{bosons, bosons}];
  (antiField[#[[1]]] = #[[2]]) & /@ fermions;
  (antiField[#[[2]]] = #[[1]]) & /@ fermions;
  (antiField[#[[1]]] = #[[2]]) & /@ complexFields;
  (antiField[#[[2]]] = #[[1]]) & /@ complexFields;
  
  ]
(* default error handler *)

setFields[a___] := Message[setFields::syntax, a];


(* countTerms cleared for DoFun3 *)
(* count the number of terms/graphs in an expression *)

countTerms[a_op,opts___]:=1;

countTerms[a_?NumericQ,___]:=0;

countTerms[a_Plus|a_Times,opts___]:=Count[a,op[___],Infinity];

countTerms[a___]:=Message[countTerms::syntax,a];



(* types of diagrams and how to extract them: connected/disconnected, 1PI/non1PI *)

(* Determines if graph is connected. *)
connectedQ[a_?NumericQ exp_op]:=connectedQ[exp]
connectedQ[graph_op] := Module[{extLegs, graphDistance},
  
  (* get external legs *)
  extLegs = 
   Cases[vertexDummies[graph], {__, a_String} /; 
     StringMatchQ[a, " leg*"], Infinity];
  
  (* take all possible distinct combinations of external legs *)
  extLegs = 
   DeleteCases[
    Union[Flatten[Outer[List, extLegs, extLegs, 1], 1]], {a_, a_}, 
    2];
    
  (* get GraphDistance for all combinations and take the maximum, 
  since infinity corresponds to disconnected *)
  graphDistance = 
   Max[GraphDistance[
       vertexDummies[graph][[All, 1]] /. Rule :> UndirectedEdge, #[[
        1]], #[[2]]] & /@ extLegs];
       
  (* return if connected or not *)
  If[graphDistance == Infinity, False, True]
  
]


(* Determine if graph is disconnected. *)
disconnectedQ[a_?NumericQ exp_op]:=disconnectedQ[exp]
disconnectedQ[graph_op] := Not[connectedQ[graph]]


(* Extract all connected diagrams. *)
getConnected[a_?NumericQ, exp_op] := a getConnected[exp]
getConnected[exp_op] /; disconnectedQ[exp] := 0
getConnected[exp_op] := exp
getConnected[exp_Plus | exp_List] := Select[exp, connectedQ]


(* Extract all disconnected diagrams. *)
getDisconnected[a_?NumericQ, exp_op] := a getDisconnected[exp]
getDisconnected[exp_op] /; connectedQ[exp] := 0
getDisconnected[exp_op] := exp
getDisconnected[exp_Plus | exp_List] := Select[exp, disconnectedQ]


(* Determine if graph is 1PI. *)
onePIQ[a_?NumericQ exp_op]:=onePIQ[exp]
onePIQ[exp_op] := Module[{props, tot},
	
  (* extract all propagators *)
  props = Cases[exp, P[___]];
  
  (* delete one propagator and check if it is connected *)
  
  tot = And @@ (connectedQ[DeleteCases[exp, #]] & /@ props);
  
  (* return if 1PI or not *)
  If[tot, True, False]
]


(* Extract all 1PI diagrams. *)
get1PI[a_?NumericQ exp_op] := a get1PI[exp]
get1PI[exp_op] /; onePIQ[exp] := exp
get1PI[exp_op] := 0
get1PI[exp_Plus | exp_List] := Select[exp, onePIQ]


(* Extract all non1PI diagrams. *)
getNon1PI[a_?NumericQ exp_op] := a getNon1PI[exp]
getNon1PI[exp_op] /; onePIQ[exp] := 0
getNon1PI[exp_op] := exp
getNon1PI[exp_Plus | exp_List] := Select[exp, Not@onePIQ[#] &]


(* predicates for fields *)

(* Note: Setting the Head of a field to 'field' does not work with pattern matching, thus fieldQ has to be used. *)
fieldQ[a_]:=Or@@(fieldType@a===# &/@$fieldTypes);

fermionQ[{}] := False (* no field *)

fermionQ[a_List]:=fermionQ[a[[1]]];(* field given together with index*)

fermionQ[a_]:=fieldType@a===fermion;

antiFermionQ[a_List]:=antiFermionQ[a[[1]]];(* field given together with index*)

antiFermionQ[a_]:=fieldType@a===antiFermion;

bosonQ[a_List]:=bosonQ[a[[1]]];(* field given together with index*)

bosonQ[a_]:=fieldType@a===boson;

complexFieldQ[a_List]:=complexFieldQ[a[[1]]];(* field given together with index*)

complexFieldQ[a_]:=fieldType@a===complex;

antiComplexFieldQ[a_List]:=antiComplexFieldQ[a[[1]]];(* field given together with index*)

antiComplexFieldQ[a_]:=fieldType@a===antiComplex;
  
(* set (anti-)commutating property: Done in setFields. *)

grassmannQ[a_List] := grassmannQ[a[[1]]](* field given together with index*)




(* ::Section:: *)
(* syntax information; experimental *)


(*SyntaxInformation[doRGE]={"ArgumentsPattern"->{_,_,OptionsPattern[]},"ArgumentsPattern"->{_,_,_,OptionsPattern[]}};*)




End[];

EndPackage[]
