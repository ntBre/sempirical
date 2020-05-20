import numpy as np
import sys

totstructs = int(sys.argv[1])
energies=[]
for structure in range(totstructs):
    file = format(structure+1, '05d')
    with open(f'structureFolder/Structure{file}.aux') as f:
        lines = f.readlines()
        for line in lines:
            if 'TOTAL_ENERGY' in line:
                energy = float(line.replace('TOTAL_ENERGY:EV=', '').replace('D', 'E'))
                print(energy)
#energies.append(float(energy))

#energies = np.array(energies)

#for item in energies:
#    print(item-np.min(energies))
