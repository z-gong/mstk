_atomic_number = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5,
                  'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
                  'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
                  'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
                  'Ti': 22, 'Fe': 26, 'Zn': 30, 'Se': 34, 'Br': 35,
                  'Kr': 36, 'Mo': 42, 'Ru': 44, 'Sn': 50, 'Te': 52,
                  'I': 53, 'Xe': 54, 'UNK': -1, 'DP': -2, 'VS': -3}

_atomic_mass = {'H': 1.008, 'He': 4.003, 'Li': 6.941, 'Be': 9.012, 'B': 10.811,
                'C': 12.011, 'N': 14.006, 'O': 15.999, 'F': 18.998, 'Ne': 20.180,
                'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.086, 'P': 30.974,
                'S': 32.065, 'Cl': 35.453, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078,
                'Ti': 47.867, 'Fe': 55.845, 'Zn': 65.38, 'Se': 78.971, 'Br': 79.904,
                'Kr': 83.798, 'Mo': 95.96, 'Ru': 101.07, 'Sn': 118.710, 'Te': 127.60,
                'I': 126.904, 'Xe': 131.293, 'UNK': 0.0, 'DP': 0.0, 'VS': 0.0}

_atomic_symbol = {v: k for k, v in _atomic_number.items()}


class Element():
    def __init__(self, arg):
        if isinstance(arg, int):
            self.number = arg
            self.symbol = _atomic_symbol[arg]
        elif isinstance(arg, str):
            self.symbol = arg
            self.number = _atomic_number[arg]
        else:
            raise Exception('Element should be initiated with atomic number or symbol')

        self.mass = _atomic_mass[self.symbol]

    def __repr__(self):
        return f'<Element: {self.symbol}>'

    @staticmethod
    def guess_from_atom_type(type):
        '''
        Guess the element from the first two characters of atom type
        Will only guess TEAM and CL&P atom types
        TEAM atom types are lowercase
        CL&P atom types are title case
        '''
        if type[0].islower():
            type = str.upper(type[0]) + type[1:]
        if type[:2] in _atomic_number.keys():
            symbol = type[:2]
        elif type[0] in _atomic_number.keys():
            symbol = type[0]
        else:
            symbol = 'UNK'
        return Element(symbol)
