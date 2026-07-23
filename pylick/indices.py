import numpy as np

from pylick._config import __table__


class IndexLibrary():
    """ Collection of Spectral indices """
    def __init__(self, file_table='./aux/tableall.dat', key_list=None):
        self.source = file_table
        self._data  = self._read_index_list(file_table)

        if key_list is not None:
            self._data  = self.__getitem__(key_list)

        self.regions    = np.array([[v[j] for j in ['blue', 'centr', 'red']] for k, v in self._data.items()])
        self.names      = np.array([v[j] for j in ['name'] for k, v in self._data.items()])
        self.units      = np.array([v[j] for j in ['unit'] for k, v in self._data.items()])
        self.tex_names  = np.array([v[j] for j in ['tex_name'] for k, v in self._data.items()])
        self.tex_units  = np.array([v[j] for j in ['tex_unit'] for k, v in self._data.items()])


    @classmethod
    def _read_index_list(self, file_table='./aux/tableall.dat'):
        """ Read index table
        Parameters
        ----------
        file_table: str
            path to the file containing the indices' definitions

        Returns
        -------
        data: dict
            dictionary of indices
            name: (name, unit, centr, band, blue, red, tex_name, tex_unit)
        """

        data = {}
	#file_table = './libs/table30.dat'
        with open('./aux/tableall.dat', 'r') as f:
            data = {}
            for line in f:
                if line[0] != '#':
                    l = line.split()
                    attr = dict(
                        centr       = (float(l[3]), float(l[4])),
                        blue        = (float(l[5]), float(l[6])),
                        red         = (float(l[7]), float(l[8])),
                        name        = l[1],
                        unit        = l[2],
                        tex_name    = l[9],
                        tex_unit    = l[10]
                    )
                    ID = int(l[0])
                    data[ID] = attr
        return data

    def __repr__(self):
        return "IndexLibrary(file_table='{:s}')".format(self.source)

    def __str__(self):
        """ Display the full library """
        msg = "{:4s}{:9s}{:10s}{:19s}{:19s}{:19s}{:10s}\n".format(
            "ID", "name", "unit", "blue", "centr", "red", "tex_name")
        msg = msg + '-'*len(msg) + '\n'
        for key, val in self._data.items():
            msg = msg + "{:<4d}{:9s}{:10s}{:8.3f}-{:8.3f}  {:8.3f}-{:8.3f}  {:8.3f}-{:8.3f}  {:s}\n".format(
                key,val['name'],val['unit'],val['blue'][0],val['blue'][1],val['centr'][0],val['centr'][1],
                val['red'][0],val['red'][1],val['tex_name'])
        return msg  

    def __len__(self):
        """ Size of the library """
        return len(self._data)

    def __getitem__(self, key_list):
        d = self._data
        return {k:d[k] for k in key_list if k in d}

    @property
    def show(self):
        return self.__str__()

    def find(self, name):
        r = []
        for item in self._data.items():
            tabname = item[1]['name']
            if name in tabname:
                r.append(tabname)
        return r
