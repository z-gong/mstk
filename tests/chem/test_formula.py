from mstools.chem.formula import Formula


def test_formula():
    formula = Formula('C10(H2C1Cl3)2')
    assert formula.atoms == {'C' : 12,
                                'H' : 4,
                                'Cl': 6
                                }
    assert formula.to_str() == 'C12H4Cl6'

    formula = Formula('UukA2(C2(H1Cl3H1)2)1O')
    assert formula.atoms == {'Uuk': 1,
                             'A'  : 2,
                             'C'  : 2,
                             'H'  : 4,
                             'Cl' : 6,
                             'O'  : 1
                             }
    assert formula.to_str() == 'C2H4A2Cl6OUuk'
