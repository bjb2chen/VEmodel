

import numpy as np


from vmio.vibronic import vIO, VMK


A = 15
N = 6
order = 2

shape = vIO.soc_model_shape_dict(A, N)
print("Example shape\n")
for k, v in shape.items(): print(k,v)
print('\n')


example_soc_model = vIO.soc_model_zeros_template_json_dict(A, N, highest_order=order)


# print(f"\nExample {order=} model\n")
# for k, v in example_soc_model.items():
#     if isinstance(v, np.ndarray):
#         print(k, v.shape)
#     else:
#         print(k, v)

# fill with necessary values
# example_soc_model[VMK.w].fill(0.3)
# example_soc_model[VMK.E].fill(0.1)
# example_soc_model[VMK.etdm].fill(0.1)
# example_soc_model[VMK.mtdm].fill(0.1)

for k in VMK.soc_coupling_list():
    if k in example_soc_model:
        example_soc_model[k].fill(3+1j)

vIO.save_model_to_JSON("./example_soc_model.json", example_soc_model)

loaded_model = vIO.load_model_from_JSON("./example_soc_model.json")

print("-"*60)
vIO.print_model_compact(loaded_model, highest_order=order)
