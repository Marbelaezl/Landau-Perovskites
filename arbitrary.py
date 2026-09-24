
# arbitrary.py - Generate solution for arbitrary quoefficients 
import numpy as np
import functions

# Initial configurarion
#               a+   T+   b+   a-   T-   b-    c    k  la+   la-  
vec = np.array([ 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 1.00000000e+00,
 8.78169796e+02, 2.51414754e+03, 0.00000000e+00, 1.75065763e+03,
 0.00000000e+00, 5.54028117e+02, 0.00000000e+00, 4.19919067e+01,
 0.00000000e+00
]) # kp



results = [[0, 0, 0, 0,0]]
status = [0]
params, couplings = functions.GenParams(vec,model="2t66p")
index = 0

#  CICLO DE SIMULACIÓN
# np.linspace(inicio, fin, pasos)
for i in np.linspace(300, 900, 1201):
    print("T = ", i)
    # Define variables: [q+, q-, T, epsilon]
    variables = np.array([
        0,
        #np.max([np.random.rand()*2, results[index][0]]),
        np.max([np.random.rand()*2, results[index][1]]),
        i,
        np.max([np.random.rand()*2, results[index][3]]),
        0
    ])

    # Energy minimization
    prov = functions.MinimizeEnergy(params, variables, couplings, mask=np.array([1, 1, 0, 1]), verbose=False, cycles=5000, delta=1e-9)

    results.append(prov[1])
    status.append(prov[0])
    index += 1

#DATA PROCESSING AND STORAGE
results_np = np.array(results)
status_np = np.array(status).reshape(-1, 1) # reshape into column

# Calculate free energy per iteration
energy_np = np.zeros((len(results_np), 1))
for i in range(len(results_np)):
    energy_np[i] = functions.FreeEnergy(params, results_np[i], couplings)


# Col 0:q+ | Col 1:q- | Col 2:T | Col 3:epsilon | Col 4:status | Col 5:free_energy
datos_totales = np.hstack([results_np, status_np, energy_np])

# GENERATE HEADER
# Build a string containing parameters in vec
header_txt = (
    f"Simulation parameters:\n"
    f"# a+={vec[0]}, T+={vec[1]}, b+={vec[2]}, a-={vec[3]}, T-={vec[4]}\n"
    f"# b-={vec[5]}, c={vec[6]}, k={vec[7]}, lambda+={vec[8]}, lambda-={vec[9]}\n"
    f"# Columnas: q+  q-  T  epsilon  status  free_energy"
)

# Save with 8 decimals (%2.8f)
np.savetxt('Sr2CuWO6-hires-biased.txt', datos_totales, header=header_txt, fmt='%2.8f')

print("--------------------------------------------------")
print("arbitraty.py finished succesfully.")
print("--------------------------------------------------")