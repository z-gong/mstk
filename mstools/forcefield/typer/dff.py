from .typer import Typer


class DffDefaultTyper(Typer):
    def __init__(self):
        super().__init__()

    def type_molecule(self, mol):
        pass


class DffDefTyper(Typer):
    def __init__(self, file):
        super().__init__()
        self._file = open(file)
        self.parse_tree()
        self._file.close()

    def parse_tree(self):
        pass


class DffExtTyper(Typer):
    def __init__(self, file):
        super().__init__()
        self._file = open(file)
        self.parse_tree()
        self._file.close()

    def parse_tree(self):
        pass
