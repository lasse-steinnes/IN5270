### Use the FEM_P2_solver class here ###
from FEM_P2 import FEM_P2_solver

Ne = 4
C = 1
D = 1

derive_b = FEM_P2_solver(C,D,Ne)
#derive_b.find_b_ref()
#derive_b.setup_element_matrix()
#derive_b.assemble()
