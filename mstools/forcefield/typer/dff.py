from .typer import Typer


class DffDefaultTyper(Typer):
    def __init__(self):
        super().__init__()


class DffDefTyper(Typer):
    def __init__(self, file):
        super().__init__()
        self._file = file
        self._parse(file)

    def _parse(self, file):
        pass


class DffExtTyper(Typer):
    def __init__(self, file):
        super().__init__()
        self._file = file
        self._parse(file)

    def _parse(self, file):
        pass
