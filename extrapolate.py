import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

#Read ionic radii information from Shannon (1976), Acta Cryst. A32, 751.
A_sites ={}
B_sites={}
def build_entry(line):
    spaces=0
    key=""
    site=""
    val=""
    for i in line:
        if i == " ":
            spaces+=1
        else:
            if spaces ==0:
                key=key+i
            if spaces==1:
                key= i+key
            if spaces==2:
                site=site+i
            if spaces==3:
                val=val+i
    if site == "12":
        A_sites[key] = val
    elif site=="6":
        B_sites[key] = val
Osize= 1.21 #Cation handled as special case since all analyzed perovskites are A2BB'O6              
        
with open("ionradii.txt", "r") as f:
    for line in f:
        build_entry(line.strip())

#Read experimental data
data_str=np.genfromtxt("extrapolation.txt",delimiter=",",dtype=str)
data_str=np.strings.strip(data_str)
print(data_str)
print(A_sites)
data_geom=[]
for i in data_str:
    prov=[]
    prov.append(A_sites[i[0]])
    prov.append(B_sites[i[1]])
    prov.append(B_sites[i[2]])
    prov.append(float(i[4]))
    data_geom.append(prov)
data_geom=np.array(data_geom,dtype=float)
print(data_geom)
t= ( (data_geom[:,0] + Osize))/(((data_geom[:,1]+data_geom[:,2])/2 + Osize)*np.sqrt(2))
print(t)
fig,ax=plt.subplots()
ax.scatter(t,data_geom[:,3])

data_trans=np.vstack([t,data_geom[:,3]]).T
print(data_trans)
import sklearn as sk

#Generate synthetic data by randomly assigning continuous transitions as either second order or tricritical
data_certain=data_geom[np.where(data_geom[:,3]!=2.5)]
print(len(data_certain))
fig,ax=plt.subplots()
data_uncertain=data_geom[np.where(data_geom[:,3]==2.5)]
#instead of guessing 50/50 split, estimate population probability of tricritic as sample probability
p=float(len(data_geom[np.where(data_geom[:,3]==3)]))/len(data_certain[np.where(data_certain[:,3] !=1)])
tries=[]
overfit_check=[]
most_common=[]
full=False

ntries=1
mpl.rcParams.update({'font.size': 22})
fig,ax=plt.subplots()
ax.hist(t,edgecolor="black",facecolor="blue", linewidth=2,rwidth=0.9,bins=15,weights=np.ones_like(t)/len(t),range=(0.9,1.05))
ax.set_xlabel("t")
ax.set_ylabel("fraction")
for i in range(0,ntries):
    rng = np.random.default_rng()
    extra=rng.binomial(1,p,size=len(data_uncertain))
    sim_order = 2.5 + 0.5*extra
    data_uncertain[:,3] = sim_order
    data_simulated= np.vstack([data_certain,data_uncertain])
    data_simulated=np.vstack([data_certain])
    #In a tricritical case there is a first and a second order transition, the probability of selecting the second order transition is 50%
    
    data_simulated[:,3]-=1
    data_simulated[:,3] *= 2
    data_simulated[np.where(data_simulated[:,3]==4),3]= 1
    #clf=sk.linear_model.LogisticRegression(C=np.inf,max_iter=1000)
    #clf=sk.naive_bayes.CategoricalNB()
    #clf=sk.ensemble.GradientBoostingClassifier(max_depth=3)
    clf=sk.ensemble.RandomForestClassifier(max_depth=3)
    #clf=sk.dummy.DummyClassifier()
    t= ( (data_simulated[:,0] + Osize))/(((data_simulated[:,1]+data_simulated[:,2])/2 + Osize)*np.sqrt(2))
    X, Y = sk.utils.resample(np.array([t]).reshape(-1,1),data_simulated[:,3])
    #X = np.array([t]).reshape(-1, 1)
    #Y = data_simulated[:,3].astype(int)

    clf.fit(X,Y)    
    tries.append(clf.score(X,Y))
    arr= np.array([np.count_nonzero(data_simulated[:,3]==0),
          np.count_nonzero(data_simulated[:,3]==1),
          np.count_nonzero(data_simulated[:,3]==2)],dtype=float)
    
    extra=rng.binomial(1,p,size=len(data_uncertain))
    sim_order = 2.5 + 0.5*extra
    data_uncertain[:,3] = sim_order
    data_simulated= np.vstack([data_certain,data_uncertain])
    data_simulated=np.vstack([data_certain])
    Y = data_simulated[:,3].astype(int)
    overfit_check.append(clf.score(X,Y))
    if (int(100*i/ntries)  < int(100*(i+1)/ntries)) and (int(100*(i+1)/ntries) % 5 ==0):
        print(str(int(100*(i+1)/ntries)) + "% complete")

tries=np.array(tries)
most_common=np.array(most_common)
overfit_check=np.array(overfit_check)
print(np.average(tries))
print(np.std(tries))
print(np.average(overfit_check))
print(np.std(overfit_check))



for i in [1.0,2.0,2.5,3.0]:
    print(len(np.where(data_trans[:,1]==i)[0]))

disc=data_trans[np.where(np.any([data_trans[:,1]==1],axis=0))]
not_disc = data_trans[np.where(data_trans[:,1]==2)]
disc[:,1]=1
not_disc[:,1]=0
not_disc = not_disc[~np.isin(not_disc[:,0], disc[:,0])]
data_disc=np.vstack([disc,not_disc])

cont= data_trans[np.where(np.any([data_trans[:,1]==2],axis=0))]
not_cont = data_trans[np.where(data_trans[:,1]==1)]
#not_cont = data_trans[np.where(data_trans[:,1]==2)]
cont[:,1]=1
not_cont[:,1]=0
not_cont = not_cont[~np.isin(not_cont[:,0], cont[:,0])]
data_cont=np.vstack([cont,not_cont])

tri= data_trans[np.where(data_trans[:,1]==3)]
not_tri = data_trans[np.where(data_trans[:,1]!=3)]
#CHECK THIS LINE
#not_tri = data_trans[np.where(data_trans[:,1]==1)]
tri[:,1]=1
not_tri[:,1]=0
not_tri = not_tri[~np.isin(not_tri[:,0], tri[:,0])]
data_tri=np.vstack([tri,not_tri])

fig,ax=plt.subplots()
ax.scatter(disc[:,0],np.ones_like(disc[:,0]))
ax.scatter(not_disc[:,0],np.zeros_like(not_disc[:,0]))
ax.set_title("disc")

fig,ax=plt.subplots()
ax.scatter(cont[:,0],np.ones_like(cont[:,0]))
ax.scatter(not_cont[:,0],np.zeros_like(not_cont[:,0]))
ax.set_title("cont")

fig,ax=plt.subplots()
ax.scatter(tri[:,0],np.ones_like(tri[:,0]))
ax.scatter(not_tri[:,0],np.zeros_like(not_tri[:,0]))
ax.set_title("tri")



# print("binary: ")
# #clf=sk.linear_model.LogisticRegression(C=np.inf,max_iter=1000)
# #clf=sk.naive_bayes.BernoulliNB()
# clf=sk.ensemble.GradientBoostingClassifier(max_depth=3)
# #clf=sk.ensemble.RandomForestClassifier(max_depth=3)
# #clf=sk.dummy.DummyClassifier()
# tries=[]
# overfit_check=[]
# for i in range(0,ntries):
#     if (int(100*i/ntries)  < int(100*(i+1)/ntries)) and (int(100*(i+1)/ntries) % 5 ==0):
#         print(str(int(100*(i+1)/ntries)) + "% complete")
#     x,y=sk.utils.resample(data_cont[:,[0]], data_cont[:,1])
#     clf.fit(x,y)
#     tries.append(clf.score(x,y))
#     x,y=sk.utils.resample(data_cont[:,[0]], data_cont[:,1])
#     overfit_check.append(clf.score(x,y))
# tries=np.array(tries)
# overfit_check=np.array(overfit_check)
# print(np.average(tries))
# print(np.std(tries))
# print(np.average(overfit_check))
# print(np.std(overfit_check))
    
