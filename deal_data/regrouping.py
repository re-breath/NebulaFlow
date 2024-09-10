import numpy as np
import ase.io as ase
import sys

def process_args():

    if len(sys.argv) < 3:
        raise ValueError("At least three arguments are required: xyzfile, direction, and at least one weight.")

    xyzfile = sys.argv[1]
    direction = sys.argv[2]
    valid_directions = ['x', 'y', 'z']
    if direction not in valid_directions:
        raise ValueError(f"Invalid direction '{direction}'. Direction must be one of {valid_directions}.")

    split_weight = [float(weight) for weight in sys.argv[3:]]
    split_num = sum(split_weight)
    return xyzfile, direction, split_weight, split_num


xyzfile, direction, split_weight, split_num = process_args()

print("\n------->initializing")
print(f"Reading : {xyzfile}...\ngrouping by : {direction}")

# # test
# # xyzfile = 'POSCAR'
# # direction =  'z'
# # split_weight = [1,2,3,1]
# # split_num = sum(split_weight)


xyzinfo = ase.read(xyzfile,index=0)


cell_x = xyzinfo.cell[0][0]
cell_y = xyzinfo.cell[1][1]
cell_z = xyzinfo.cell[2][2]
positions = xyzinfo.get_positions()
split_weight = np.array(split_weight)

match direction:
    case 'z':
        pos = positions[:, 2]
        cell_length = cell_z
    case 'y':
        pos = positions[:, 1]
        cell_length = cell_y
    case 'x':
        pos = positions[:, 0]
        cell_length = cell_x
    case _:
        print("Invalid direction")
        raise ValueError("Invalid direction specified.")

split_zone = np.cumsum(split_weight)
bins = cell_length/(split_num)*split_zone
bins = bins[:-1]
print("group flag : ",bins)

indices = np.digitize(pos,bins)
xyzinfo.arrays['group'] = indices

print("\n--------->ending")
group_id,group_count = np.unique(indices,return_counts=True)
for i in group_id:
    print(f"Group {i}: {group_count[i]}")


ase.write('grouped.xyz',xyzinfo)
print("\nGroup information added successfully.\n")