class Formula:
    def __init__(self, formula: str):
        chars = self.extract_chars(formula)
        chars = self.expand_chars(chars)
        self.atoms = self.count_atoms(chars)

    @staticmethod
    def extract_chars(formula):
        chars = []
        i = 0
        while i < len(formula):
            c = formula[i]
            if c.isupper():
                chars.append(c)
                i += 1
                while i < len(formula):
                    if formula[i].islower():
                        chars[-1] += formula[i]
                        i += 1
                    else:
                        break
            elif c.isdigit():
                chars.append(c)
                i += 1
                while i < len(formula):
                    if formula[i].isdigit():
                        chars[-1] += formula[i]
                        i += 1
                    else:
                        break
            elif c == '(' or c == ')':
                chars.append(c)
                i += 1
            else:
                raise Exception('Invalid character: %s' % c)

        if chars[0].isdigit():
            raise Exception('Invalid formula')
        if chars.count('(') != chars.count(')'):
            raise Exception('Unmatched brackets')
        return chars

    @staticmethod
    def expand_chars(chars):
        chars = chars[:]
        temp_expanded_chars = []
        temp_chars = []
        while '(' in chars:
            _start = False
            i = 0
            while i < len(chars):
                c = chars[i]
                if c == '(':
                    _start = True
                    temp_expanded_chars = chars[:i]
                    temp_chars = []
                    i += 1
                elif c == ')':
                    if not _start:
                        raise Exception('Unmatched brackets')
                    if i + 1 < len(chars) and chars[i + 1].isdigit():
                        n = int(chars[i + 1])
                        temp_expanded_chars += temp_chars * n
                        if i + 2 < len(chars):
                            temp_expanded_chars += chars[i + 2:]
                    else:
                        n = 1
                        temp_expanded_chars += temp_chars * n
                        if i + 1 < len(chars):
                            temp_expanded_chars += chars[i + 1:]
                    chars = temp_expanded_chars
                    break
                else:
                    if _start:
                        temp_chars.append(c)
                    i += 1

        return chars

    @staticmethod
    def count_atoms(chars):
        counts = {}
        for i in range(len(chars)):
            c = chars[i]
            if c.isdigit():
                counts[chars[i - 1]] += int(c) - 1
            else:
                if c not in counts:
                    counts[c] = 1
                else:
                    counts[c] += 1
        return counts


s = 'UukA2(C2(H1Cl3H1)2)1O'
formula = Formula(s)
print(formula.atoms)
