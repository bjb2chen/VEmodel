import os

filename = "auto"
#filename = "auto_h2o_FC_linear_PBF100_tf50.00"
#filename = "auto_h2o_vibronic_linear_PBF100_tf50.00"
#filename = "auto_h2o_vibronic_quadratic_PBF100_tf50.00"
output_file = f"{filename}"
input_file = f"{filename}"
density = 4000
left_eV = 11
right_eV = 21
tau = 40
iexp = 1

command = (
	"autospec84 "
	f"-o {output_file:s} "
	f"-f {input_file:s} "
	f"-p {density:d} "
	f"{left_eV} {right_eV} eV "
	f"{tau:d} "
	f"{iexp:d} "
	)

print(command)
os.system(command)
