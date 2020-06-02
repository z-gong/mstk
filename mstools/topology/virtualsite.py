from .atom import Atom


def create_virtual_site(type: str, parents: [Atom], parameters: [float]):
    if type == TwoLineSite.__name__:
        return TwoLineSite(parents, parameters)
    elif type == ThreePlaneSite.__name__:
        return ThreePlaneSite(parents, parameters)
    else:
        raise Exception('Unknown virtual site type')


class VirtualSite():
    def __init__(self, type: str, parents: [Atom], parameters: [float]):
        self.type = type
        self.parents = parents[:]
        self.parameters = parameters[:]


class TwoLineSite(VirtualSite):
    def __init__(self, parents, parameters):
        if len(parents) != 2 or len(parameters) != 1:
            raise Exception('Invalid number of parents or parameters')

        super().__init__(self.__class__.__name__, parents, parameters)

    def calc_position(self):
        if not all(atom.has_position for atom in self.parents):
            raise Exception('Position for parent atoms are required')

        a1, a2 = self.parents
        p1 = self.parameters
        return a1.position * p1 + a2.position * (1 - p1)


class ThreePlaneSite(VirtualSite):
    def __init__(self, parents, parameters):
        if len(parents) != 3 or len(parameters) != 2:
            raise Exception('Invalid number of parents or parameters')

        super().__init__(self.__class__.__name__, parents, parameters)

    def calc_position(self):
        if not all(atom.has_position for atom in self.parents):
            raise Exception('Position for parent atoms are required')

        a1, a2, a3 = self.parents
        p1, p2 = self.parameters
        return a1.position * p1 + a2.position * p2 + a3.position * (1 - p1 - p2)
