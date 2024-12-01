## This is an example for defining atom types
## It is not supposed to be used for production

TypingEngine SmartsTyper

TypeDefinition

H       [#1]
H1      [#1;X1]
H1p     [#1;X1][#7,#8]
C       [#6]
C4      [#6;X4]
C3      [#6;X3]
C3a     [c;X3]
C3=O    [C;X3]=O
C2      [#6;X2]
N       [#7]
N3      [#7;X3]
N3CO    [N;X3]C=O
O       [#8]
O2      [#8;X2]
O2w     [O;X2;H2]
O1      [O;X1]


HierarchicalTree

H
    H1
        H1p
C
    C4
    C3
        C3a
        C3=O
    C2
N
    N3
        N3CO
O
    O2
        O2w
    O1
