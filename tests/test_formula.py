from mstools.formula import Formula

def test_read():
    formula = Formula('C10(H2C1Cl3)2')
    print(formula.atomdict)