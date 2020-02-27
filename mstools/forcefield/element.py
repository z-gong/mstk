atomic_number = {'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5,
                 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
                 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
                 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
                 'Ti': 22, 'Fe': 26, 'Zn': 30, 'Se': 34, 'Br': 35,
                 'Kr': 36, 'Mo': 42, 'Ru': 44, 'Sn': 50, 'Te': 52,
                 'I': 53, 'Xe': 54, 'UNK': -1, 'D': -2}

atomic_symbol = dict([(v, k) for k, v in atomic_number.items()])


class Element():
    def __init__(self, number):
        self.number = number
        self.symbol = atomic_symbol[number]

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
        if type[:2] in atomic_number:
            number = atomic_number[type[:2]]
        elif type[0] in atomic_number:
            number = atomic_number[type[0]]
        else:
            print('warning: unknown element for atom type ' + type)
            number = -1
        return Element(number)
