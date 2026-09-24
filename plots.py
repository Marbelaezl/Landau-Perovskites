import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
# 1. Cargamos el archivo único de 6 columnas
data = np.loadtxt('Sr2NiMoO6-hires-biased.txt')
data2 = np.loadtxt('Sr2NiMoO6-hires-unbiased.txt')

# data3 = np.loadtxt('La2CoMnO6-hires-biased.txt')
# data4 = np.loadtxt('La2CoMnO6-hires-unbiased2.txt')
# data[np.where(data[:,2] < 450)]=data3[np.where(data[:,2] < 450)]
#data2[np.where(data2[:,2] < 450)]=data4[np.where(data2[:,2] < 450)]
#data2 = np.loadtxt('datos_simulacion2.txt')
os.chdir('./exp_data/processed/Sr2NiMoO6/')
merged=np.genfromtxt("Sr2NiMoO6-merged.txt")
merged2=np.genfromtxt("Sr2NiMoO6-merged-certain.txt")

mpl.rcParams.update({'font.size': 16})
prov= merged[np.where(merged[:,3]==0)]

fig,ax=plt.subplots()
# ax.plot(data2[:,2],data2[:,6])
# ax.plot(data[:,2],data[:,6])
# ax.plot(data3[:,2],data3[:,6])
# ax.plot(data4[:,2],data4[:,6])

#vref=150.918123

#vref=np.max(merged[:,1])
print(merged)
vref=123.224023
merged[:,1] = (merged[:,1]/vref) - 1
merged2[:,1] = (merged2[:,1]/vref) - 1

# 2. SEPARACIÓN DE DATOS (Mapeo de las 6 columnas)
# Col 0:q+, Col 1:q-, Col 2:T, Col 3:epsilon, Col 4:status, Col 5:energy
results = data[:, 0:4]  # Toma las primeras 4 columnas (q+, q-, T, epsilon)
results2 = data2[1:, 0:4]  # Toma las primeras 4 columnas (q+, q-, T, epsilon)
status  = data[:, 4]    # Toma la quinta columna (status de error)
energy  = data[:, 5]    # Toma la sexta columna (energía libre)

results=results[1:]
#results[:,0]=0
#results2[:,0]=0
vec = ([-1, 1, 0.3, -1.5, 5, 0.6, -1.0, 0.25, 0.1, -0.2])

# --- BLOQUE DE GRÁFICAS ORIGINAL ---
fig, ax = plt.subplot_mosaic([[0, 2],
                               [0, 2]],
                              figsize=(9, 6), layout="constrained")

# Identificamos errores usando la columna status que ya cargamos
errors = results[np.where(status != 0)]
results[:,[0,1]] = np.abs(results[:,[0,1]])
results2[:,[0,1]] = np.abs(results2[:,[0,1]])
# Gráfica q+ (Columna 2 de results es T, Columna 0 es q+)
# ax[1].plot(results[:,2], results[:,0],color="blue",label="Full-data estimation")
# ax[1].plot(results2[:,2], results2[:,0],color="red",label="Reduced bias estimation",linestyle="--")
# #ax[1].scatter(errors[:,2], errors[:,0], color="red",s=4)
# ax[1].set_ylabel(r'$q_+ (\AA)$')

#Gráfica q- (Columna 1 de results es q-)
ax[2].plot(results[:,2], results[:,1],color="blue",label="Full-data estimation")
ax[2].plot(results2[:,2], results2[:,1],color="red",label="Reduced bias estimation",linestyle="--")
#ax[2].scatter(errors[:,2], errors[:,1], color="red",s=4)
ax[2].set_ylabel(r'$q_- (\AA)$')

# Gráfica epsilon (Columna 3 de results es epsilon)
ax[0].plot(results[:,2], results[:,3],color="blue",label="Full-data estimation")
ax[0].plot(results2[:,2], results2[:,3],color="red",label="Reduced bias estimation",linestyle="--")
#ax[0].scatter(errors[:,2], errors[:,3], color="red",s=4)
ax[0].set_ylabel(r'$\epsilon$')

for i in range(0,3):
    if i!=1:
        ax[i].set_xlabel("T(K)")
#data=data[:,[2,3,0,1]]

for i in range(0,3):
    if i!=1:
        ax[i].scatter(merged[:,0],merged[:,i+1],color="grey",s=10)
        ax[i].scatter(merged2[:,0],merged2[:,i+1],color="black",label="data",s=10)

ax[0].legend()
#Reference vertical lines (I2/m phase)
# ax[0].plot([873, 873], [np.min(results[:,3]), np.max(results[:,3])], linestyle="--", color="blue")
# ax[0].plot([973, 973], [np.min(results[:,3]), np.max(results[:,3])], linestyle="--", color="blue")

# ax[1].plot([873, 873], [np.min(results[:,0]), np.max(results[:,0])], linestyle="--", color="blue")
# ax[1].plot([973, 973], [np.min(results[:,0]), np.max(results[:,0])], linestyle="--", color="blue")

# ax[2].plot([873, 873], [np.min(results[:,1]), np.max(results[:,1])], linestyle="--", color="blue")
# ax[2].plot([973, 973], [np.min(results[:,1]), np.max(results[:,1])], linestyle="--", color="blue")



# --- GRÁFICO DE ENERGÍA (AHORA MÁS RÁPIDO) ---
# Ya no usamos el bucle for ni functions.FreeEnergy porque ya tenemos 'energy'
# fig, ax_en = plt.subplots()
# ax_en.scatter(results[:, 2][np.where(status==0)], energy[np.where(status==0)],s=4,color="black") # results[:, 2] es la temperatura T
# ax_en.set_xlabel("T")
# ax_en.set_ylabel("Free Energy")

plt.show()