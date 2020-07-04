from .atom import Atom


class VirtualSite():
    '''
    Base class for virtual site definitions
    '''

    _class_map = {}

    def __init__(self):
        self.parents = []
        self.parameters = []

    @staticmethod
    def register(klass):
        VirtualSite._class_map[klass.__name__] = klass

    @staticmethod
    def create(type, parents, parameters):
        '''
        Factory function for creating a new virtual site

        Parameters
        ----------
        type: str
            The type of virtual site to be created
        parents : list of Atom
            The parent atoms of the new virtual site
        parameters : list of float
            Parameters for calculating the position of virtual site from parent atoms

        Returns
        -------
        virtual_site : subclass of VirtualSite
        '''

        try:
            cls = VirtualSite._class_map[type]
        except:
            raise Exception('Unknown virtual site type. Valid types: '
                            + str(list(VirtualSite._class_map.keys())))

        return cls(parents, parameters)

    def calc_position(self):
        '''
        Calculate the position of virtual site from parent atoms.

        Returns
        -------
        position : array_like
        '''
        raise NotImplementedError('This method should be implemented by subclasses')


class TwoLineSite(VirtualSite):
    '''
    A virtual site in the same line defined by two atoms

    Parameters
    ----------
    parents : list of Atom
    parameters : list of float
    '''
    def __init__(self, parents, parameters):
        if len(parents) != 2 or len(parameters) != 1:
            raise Exception('Invalid number of parents or parameters')

        super().__init__()
        self.parents = parents[:]
        self.parameters = parameters[:]

    def calc_position(self):
        if not all(atom.has_position for atom in self.parents):
            raise Exception('Position for parent atoms are required')

        a1, a2 = self.parents
        p1 = self.parameters
        return a1.position * p1 + a2.position * (1 - p1)


VirtualSite.register(TwoLineSite)


class ThreePlaneSite(VirtualSite):
    '''
    A virtual site in the same plane defined by three atoms

    Parameters
    ----------
    parents : list of Atom
    parameters : list of float
    '''
    def __init__(self, parents, parameters):
        if len(parents) != 3 or len(parameters) != 2:
            raise Exception('Invalid number of parents or parameters')

        super().__init__()
        self.parents = parents[:]
        self.parameters = parameters[:]

    def calc_position(self):
        if not all(atom.has_position for atom in self.parents):
            raise Exception('Position for parent atoms are required')

        a1, a2, a3 = self.parents
        p1, p2 = self.parameters
        return a1.position * p1 + a2.position * p2 + a3.position * (1 - p1 - p2)


VirtualSite.register(ThreePlaneSite)
