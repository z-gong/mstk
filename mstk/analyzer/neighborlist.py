import itertools
import numpy as np


class NeighborList:
    '''
    Parameters
    ----------
    box : np.ndarray of shape (3,)
    cutoff : float

    Attributes
    ----------
    box : np.ndarray of shape (3,)
    n_cell : np.ndarray of shape (3,)
    cell_size : np.ndarray of shape (3,)
    cell_indexes : list of tuple of int
    '''

    def __init__(self, box, cutoff):
        if any(box / 3 < cutoff):
            raise Exception('Cutoff larger than box/3')
        n = np.zeros(3, dtype=int)
        for i in range(3):
            while True:
                if box[i] / (n[i] + 1) < cutoff:
                    break
                n[i] += 1

        self.box = box
        self.n_cell = n
        self.cell_indexes = list(itertools.product(range(n[0]), range(n[1]), range(n[2])))
        self.cell_size = box / n
        self._cells = [[[[] for _ in range(n[2])] for _ in range(n[1])] for _ in range(n[0])]

    def build(self, positions):
        '''
        Parameters
        ----------
        positions : np.ndarray
        '''
        positions -= np.floor(positions / self.box) * self.box
        locations = np.clip(np.floor(positions / self.cell_size).astype(int), 0, self.n_cell - 1)
        for i, (x, y, z) in enumerate(locations):
            self._cells[x][y][z].append(i)

    def get_cell(self, index):
        '''
        Parameters
        ----------
        index : tuple of int

        Returns
        -------
        index_atoms : list of int
        '''
        return self._cells[index[0]][index[1]][index[2]]

    def get_interacting_cell(self, index):
        '''
        Return the 27 interacting cells as a single large cell

        Parameters
        ----------
        index : tuple of int

        Returns
        -------
        index_atoms : list of int
        '''
        idx_surround = [[], [], []]
        for i in range(3):
            idx = index[i]
            idx_surround[i] = [idx - 1 if idx - 1 >= 0 else self.n_cell[i] - 1,
                               idx,
                               idx + 1 if idx + 1 <= self.n_cell[i] - 1 else 0
                               ]
        cell = []
        for index in itertools.product(*idx_surround):
            cell.extend(self.get_cell(index))

        return cell
