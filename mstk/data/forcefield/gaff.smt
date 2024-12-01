## This is an unofficial implementation of GAFF atom types based on Table 1 in the following article:
## * Wang, J.; Wolf, R. M.; Caldwell, J. W.; Kollman, P. A.; Case, D. A. Development and Testing of a General Amber Force Field. Journal of Computational Chemistry 2004, 25 (9), 1157–1174. https://doi.org/10.1002/jcc.20035.
## It covers the 35 basic atom types and 23 special atom types.
## Because of the special treatment for conjugated atoms, it requires GaffTyper typing engine.

## RDKit is used here for determining aromaticity, which is not consistent with GAFF.
## Especially, GAFF does not consider 5-member rings as aromatic.
## For some molecules, this definition does not give the same results as the official implementation in AmberTools.
## Therefore, this definition is mainly used for optimizing a force field by taking GAFF as an initial guess.

## Hybridization assignment is nasty. E.g. N in N-c1ccccc1 is SP2 hybridized, whereas P in P-c1ccccc1 is SP3 hybridized.
## N in N-C=C is considered as SP2 hybridized by RDKit, but it is not a strong planer structure.
## Therefore, hybridization is avoided if the type can be correctly assigned with other information.

TypingEngine GaffTyper

TypeDefinition

c      [C;X3]=[O,S;X1]
c1     [C;X2]
c2     [C;X3]
c3     [C;X4]
ca     [c]
n      [N;X3][C,S,P]=[O,S]      ## amide. Cannot use hybridization. RDKit treat N in NC=O as SP2, but in NS=O as SP3
n1     [N^1;X1]                 ## SP1 nitrogen. E.g. nitrile
n2     [N^2;X2]                 ## SP2 nitrogen with 2 coordinates and a double bond. E.g. imine
n3     [N^3;X3]                 ## SP3 nitrogen with 3 coordinates. E.g. amine
n4     [N^3;X4]                 ## SP3 nitrogen with 4 coordinates. E.g. ammonium
na     [N;X3](*=*)(*=*)         ## SP2 nitrogen with 3 coordinates. E.g. amine connected two alkene groups. Cannot use N^2, because RDKit considers N in NC=C as SP2
nh     [N;X3](*)(*)[a]          ## amine connected to aromatic
no     [N;X3](=O)(~O)           ## nitro group
o      [O;X1]                   ## O in carbonyl or nitro group
oh     [O;X2;H]
os     [#8;X2;H0]               ## there is no aromatic O in GAFF
s2     [S;X1]=*
sh     [S;X2;H]
ss     [#16;X2;H0]              ## there is no aromatic S in GAFF
s4     [S;X3]=*
s6     [S;X4]=*
p2     [P;X2]=*
p3     [P;X3](*)(*)*            ## three coordinates without double bond. P in P-c1ccccc1 is SP3 hybridized
p4     [P;X3]=*
p5     [P;X4]=*
hc     [H]C
ha     [H]c
hn     [H][#7]
ho     [H][#8]
hs     [H][#16]
hp     [H][#15]
f      [F]
cl     [Cl]
br     [Br]
i      [I]
cc|cd  [C;X3](=;@*)-;@*!-;@*
ce|cf  [C;X3](=*)-*!-*
cg|ch  [C;X2](#*)-*!-*
cp|cq  [c]-[a]
cu     [C;X3;r3]
cv     [C;X3;r4]
cx     [C;X4;r3]
cy     [C;X4;r4]
nb     [n]
nc|nd  [N;X2](=;@*)-;@*!-;@*
ne|nf  [N;X2](=*)-*!-*
pb     [p]
pc|pd  [P;X2](=;@*)-;@*!-;@*
pe|pf  [P;X2](=*)-*!-*
px     [P;X3](=*)-*!-*
py     [P;X4](=*)-*!-*
sx     [S;X3](=*)-*!-*
sy     [S;X4](=*)-*!-*
h1     [H]C~[#7,#8,F,Cl,Br]
h2     [H]C(~[#7,#8,F,Cl,Br])~[#7,#8,F,Cl,Br]
h3     [H]C(~[#7,#8,F,Cl,Br])(~[#7,#8,F,Cl,Br])~[#7,#8,F,Cl,Br]
h4     [H]c~[#7,#8,F,Cl,Br]
h5     [H]c(~[#7,#8,F,Cl,Br])~[#7,#8,F,Cl,Br]


HierarchicalTree

cu
cv
cx
cy
c
c1
    cg|ch
c2
    ce|cf
        cc|cd
c3
ca
    cp|cq
n
n1
n2
    ne|nf
        nc|nd
n3
n4
no
nh
na
o
oh
os
s2
sh
ss
s4
    sx
s6
    sy
p2
    pe|pf
        pc|pd
p3
p4
    px
p5
    py
hc
    h1
        h2
            h3
ha
    h4
        h5
hn
ho
hs
hp
f
cl
br
i
nb
pb
