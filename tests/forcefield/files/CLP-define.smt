TypeDefinition
## basic
H_1 [H]
C_2 [#6X2]
C_3 [#6X3]
C_4 [#6X4]
N_1 [#7X1]
N_2 [#7X2]
N_3 [#7X3]
O_1 [#8X1]
O_2 [#8X2]
F_1 [#9X1]
S_4 [#16X4]

## alkane
HC [H][CX4]
CT [CX4;H3]
CS [CX4;H2]

## imidazolium
NA n1cncc1 1
CR c1nccn1 1
CW c1cncn1 1
HCR [H]c1nccn1 1
HCW [H]c1cncn1 1
C1 [CX4;H2,H3]n1cncc1 1
H1 [H][CX4]n1cncc1 1
C2 [CX4;H2][CX4]n1cncc1 1
CE [CX4;H3][CX4]n1cncc1 1

## benzene
CA c1ccccc1
HA [H]c1ccccc1

## ethylbenzene-imidazolium
CAT c1(C)ccccc1
CAO [cH]1c(C)cccc1
CAM [cH]1cc(C)ccc1
CAP [cH]1ccc(C)cc1
HAT [H][$(c1c(C)cccc1),$(c1cc(C)ccc1),$(c1ccc(C)cc1)]
C_4H2_CA Cc1ccccc1
HT [H][CX4]c1ccccc1
C_4H2_NA [CH2](n1cncc1)Cc1ccccc1 1

## dicyanoimide
N3A [N-]C#N -1
CZA C([N-])#N -1
NZA N#C[N-] -1

## FSI TFSI
NBT [N-](S(=O)(=O))S(=O)(=O) -1
SBT S(=O)(=O)[N-] -1
OBT O=S(=O)[N-] -1
FSI FS(=O)(=O)[N-] -1
CBT C(F)(F)(F)S(=O)(=O)[N-] -1
F1 F[CX4]

HierarchicalTree
H_1
    HCR
    HCW
    HC
        H1
        HT
    HA
        HAT
C_2
    CZA
C_3
    CR
    CW
    CA
        CAT
        CAO
        CAM
        CAP
C_4
    C1
        C_4H2_NA
    CT
        CE
    CS
        C_4H2_CA
        C2
    CBT
N_1
    NZA
N_2
    N3A
    NBT
N_3
    NA
O_1
    OBT
O_2
F_1
    FSI
    F1
S_4
    SBT
