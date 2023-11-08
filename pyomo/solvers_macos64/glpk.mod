/* Created by glpkXl, Friday, November 3, 2023 5:19:40 PM */


set materials   ;
set components   ;
param analysis {materials, components}  symbolic;

set allAtoms   ;
param massA {allAtoms}  ;

set stoechiometry   dimen 2;
set molecules   := setof{(m,a) in stoechiometry} m;
set atoms   := setof{(m,a) in stoechiometry} a;
param stoech {m in molecules, a in atoms}  default 0;
param massM {m in molecules}  :=sum {a in atoms} stoech[m,a]*massA[a];
param stoechM {m in molecules, a in atoms}  :=stoech[m,a]*massA[a]/massM[m];

set oximetry   dimen 2;
set oxides   := setof{(m,o) in oximetry} o;
param oxi {m in molecules, o in oxides}  default 0;
param oxiM {m in molecules, o in oxides}  :=oxi[m,o]*massM[o]/massM[m];

set flowsheet   dimen 3;
set from   :=setof{(f,t,p) in flowsheet} f;
set to   :=setof{(f,t,p) in flowsheet} t;
set boxes   := from inter to;
param balSign {b in boxes, (f,t,p) in flowsheet}  := if t=b then 1 else if f=b then -1 else 0;

var Q {flowsheet}  >=0;
var A {flowsheet, components}  >=0;
var O {flowsheet, oxides}  >=0;

s.t. known {(f,t,p) in flowsheet, c in components: analysis[p,c]!='x'}: Q[f,t,p] * analysis[p,c] =A[f,t,p,c];
s.t. balanceA {b in boxes, a in atoms}: sum {(f,t,p) in flowsheet, m in molecules} balSign[b,f,t,p] * A[f,t,p,m] * stoechM[m,a] =0;
s.t. oxiDef {(f,t,p) in flowsheet, o in oxides}: sum {c in molecules} A[f,t,p,c] * oxiM[c,o] =O[f,t,p,o];


s.t. dusts:  A['PH','EX','Tuffeau_67','tdry'] =0.05*A['RM','PH','Tuffeau_67','tdry'];
s.t. GJBU:  A['FU','BU','Charbon_Terval','LHV'] =2.5*A['CO4','SI','clinker','tcl'];
s.t. GJPR : A['FU','PR','Charbon_Terval','LHV'] =1.5*A['CO4','SI','clinker','tcl'];
s.t. production:  A['CO4','SI','clinker','tcl'] =100;
s.t. c2sCli:  A['CO4','SI','clinker',"c2s"] =0.20*A['CO4','SI','clinker','tcl'];
s.t. c3aCli:  A['CO4','SI','clinker',"c3a"] =0.06*A['CO4','SI','clinker','tcl'];
s.t. c4afCli:  A['CO4','SI','clinker',"c4af"] =0.09*A['CO4','SI','clinker','tcl'];

end;
