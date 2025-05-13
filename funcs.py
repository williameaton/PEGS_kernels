import ensightreader
import numpy as np
class EnData():
    def __init__(self, iproc, fpath, region=1, part=1):

        self.iproc  = iproc
        self.fpath  = fpath
        self.region = region
        self.part   = part

        self.vars    = {}
        self.vardata = {}



    def load_ProcData(self, varname):

            self.casefpath = f'{self.fpath}/reg{self.region}_proc{self.iproc}.case'
            self.case = ensightreader.read_case(self.casefpath)
            self.geofile = self.case.get_geometry_model()

            self.part_names = self.geofile.get_part_names()
            part = self.geofile.get_part_by_name(self.part_names[0])
            self.Nnodes = part.number_of_nodes

            with self.geofile.open() as fp_geo:
                self.node_coordinates = part.read_nodes(fp_geo)

            # Read data for this kernel
            variable = self.case.get_variable(varname)
            self.vars[varname] = variable
            with open(variable.file_path, "rb") as fp_var:
                self.vardata[varname] = variable.read_node_data(fp_var, 1)




def norm(x):
    return x/np.amax(np.abs(x))