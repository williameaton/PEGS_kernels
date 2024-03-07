import ensightreader as en

class ProcData():

    def __init__(self, proc, fpath, region=1):
        self.proc = proc
        self.fpath = fpath
        self.region = region

        self.data = {}
        self.case = en.read_case(f"{fpath}/reg{self.region}_proc{self.proc}.case")

        self.geofile = self.case.get_geometry_model()

        self.part_names = self.geofile.get_part_names()
        self.part = self.geofile.get_part_by_name(self.part_names[0])
        self.Nnodes = self.part.number_of_nodes

        with open(self.geofile.file_path, "rb") as fp_geo:
            self.coord = self.part.read_nodes(fp_geo)
            for block in self.part.element_blocks:
                self.connectivity = block.read_connectivity(fp_geo)


    def load_var(self, varname):
        variable =  self.case.get_variable(varname)
        with open(variable.file_path, "rb") as fp_var:
            d = variable.read_node_data(fp=fp_var, part_id=self.part.part_id)

        self.data[varname] = d

