from .typer import Typer


class ZftTyper(Typer):
    def __init__(self, file):
        super().__init__()
        self._file = open(file)
        self.parse_tree()
        self._file.close()

    def parse_tree(self):
        pass
