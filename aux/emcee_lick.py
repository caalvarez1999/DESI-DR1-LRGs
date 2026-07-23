import numpy as np
import emcee

####################################################################
#Create a dictionary linking Lick index name and identifying number
Lick_indices={"HdA":0,"HdF":1,"CN1":2,"CN2":3,"Ca4227":4,"G4300":5,"HgA":6,"HgF":7,"Fe4383":8,"Ca4455":9,"Fe4531":10,"C24668":11,"Hb":12,"Fe5015":13,"Mg1":14,"Mg2":15, "Mgb":16,"Fe5270":17,"Fe5335":18,"Fe5406":19,"Fe5709":20,"Fe5782":21, "NaD":22,"TiO1":23,"TiO2":24}
Lick_indices27={"D4000":0, "Dn4000":1, "HdA":2,"HdF":3,"CN1":4,"CN2":5,"Ca4227":6,"G4300":7,"HgA":8,"HgF":9,"Fe4383":10,"Ca4455":11,"Fe4531":12,"C24668":13,"Hb":14,"Fe5015":15,"Mg1":16,"Mg2":17, "Mgb":18,"Fe5270":19,"Fe5335":20,"Fe5406":21,"Fe5709":22,"Fe5782":23, "NaD":24,"TiO1":25,"TiO2":26}
Lick_indices29={"CaIIK":0, "CaIIH":1, "D4000":2, "Dn4000":3, "HdA":4,"HdF":5,"CN1":6,"CN2":7,"Ca4227":8,"G4300":9,"HgA":10,"HgF":11,"Fe4383":12,"Ca4455":13,"Fe4531":14,"C24668":15,"Hb":16,"Fe5015":17,"Mg1":18,"Mg2":19, "Mgb":20,"Fe5270":21,"Fe5335":22,"Fe5406":23,"Fe5709":24,"Fe5782":25, "NaD":26,"TiO1":27,"TiO2":28}
Lick_indices30={"CaIIK":0, "CaIIH":1, "D4000":2, "Dn4000":3, "HdA":4,"HdF":5,"CN1":6,"CN2":7,"Ca4227":8,"G4300":9,"HgA":10,"HgF":11,"Fe4383":12,"Ca4455":13,"Fe4531":14,"C24668":15,"Hb0":16,"Hb":17,"Fe5015":18,"Mg1":19,"Mg2":20, "Mgb":21,"Fe5270":22,"Fe5335":23,"Fe5406":24,"Fe5709":25,"Fe5782":26, "NaD":27,"TiO1":28,"TiO2":29}
Lickemission_indices={"CaIIK":0, "CaIIH":1, "D4000":2, "Dn4000":3, "HdA":4,"HdF":5,"CN1":6,"CN2":7,"Ca4227":8,"G4300":9,"HgA":10,"HgF":11,"Fe4383":12,"Ca4455":13,"Fe4531":14,"C24668":15,"Hb0":16,"Hb":17,"Fe5015":18,"Mg1":19,"Mg2":20, "Mgb":21,"Fe5270":22,"Fe5335":23,"Fe5406":24,"Fe5709":25,"Fe5782":26, "NaD":27,"TiO1":28,"TiO2":29,"OII":100,"OIII":101,"Ha":102}

Lickopen_indices={"OII":0,"CaIIK":1, "CaIIH":2, "D4000":3, "Dn4000":4, "HdA":5,"HdF":6,"CN1":7,"CN2":8,"Ca4227":9,"G4300":10,"HgA":11,"HgF":12,"Fe4383":13,"Ca4455":14,"Fe4531":15,"C24668":16,"Hb0":17,"Hb":18,"OIII5007":19,"Fe5015":20,"Mg1":21,"Mg2":22, "Mgb":23,"Fe5270":24,"Fe5335":25,"Fe5406":26,"Fe5709":27,"Fe5782":28, "NaD":29,"TiO1":30,"TiO2":31,"Ha":32}
Lickopen_type={"OII":"A","CaIIK":"A", "CaIIH":"A", "D4000":"break", "Dn4000":"break", "HdA":"A","HdF":"A","CN1":"mag","CN2":"mag",
"Ca4227":"A","G4300":"A","HgA":"A","HgF":"A","Fe4383":"A","Ca4455":"A","Fe4531":"A","C24668":"A","Hb0":"A","Hb":"A",
"OIII5007":"A","Fe5015":"A","Mg1":"mag","Mg2":"mag", "Mgb":"A","Fe5270":"A","Fe5335":"A","Fe5406":"A","Fe5709":"A","Fe5782":"A", 
"NaD":"A","TiO1":"mag","TiO2":"mag","Ha":"A"}

Lick_indicesSDSS = {"d4000":0, "d4000_n":1, "lick_hd_a":2, "lick_hd_f":3, "lick_cn1":4, "lick_cn2":5, "lick_ca4227":6, "lick_g4300":7, "lick_hg_a":8, "lick_hg_f":9, "lick_fe4383":10, "lick_ca4455":11, "lick_fe4531":12, "lick_c4668":13, "lick_hb":14, "lick_fe5015":15, "lick_mg1":16, "lick_mg2":17, "lick_mgb":18, "lick_fe5270":19, "lick_fe5335":20, "lick_fe5406":21, "lick_fe5709":22, "lick_fe5782":23, "lick_nad":24, "lick_tio1":25, "lick_tio2":26}


Lick_Knowles4 = {"Hb0":0, "Mgb":1,"Fe5270":2,"Fe5335":3,"MgFe":4, "MgFe_prima":5}


#Return Lick index name associated to a given reference number
def I_name(index_number):
    for index in Lick_indices:
        if Lick_indices[index]==index_number:
            return index

#Return Lick index name associated to a given reference number
def I_name27(index_number):
    for index in Lick_indices27:
        if Lick_indices27[index]==index_number:
            return index

#Return Lick index name associated to a given reference number
def I_name29(index_number):
    for index in Lick_indices29:
        if Lick_indices29[index]==index_number:
            return index

#Return Lick index name associated to a given reference number
def I_name30(index_number):
    for index in Lick_indices30:
        if Lick_indices30[index]==index_number:
            return index

#Return Lick index name associated to a given reference number
def I_nameall(index_number):
    for index in Lickemission_indices:
        if Lickemission_indices[index]==index_number:
            return index

#Return Lick index name associated to a given reference number
def I_nameopen(index_number):
    for index in Lickemission_indices:
        if Lickopen_indices[index]==index_number:
            return index

#Return Lick index name associated to a given reference number
def I_nameSDSS(index_number):
    for index in Lick_indicesSDSS:
        if Lick_indicesSDSS[index]==index_number:
            return index

#Return Lick index name associated to a given reference number
def I_nameKnowles4(index_number):
    for index in Lick_Knowles4:
        if Lick_Knowles4[index]==index_number:
            return index







#Return index reference number associated to a given Lick index
def I_number(index_name):
    index_numbers=[]
    for index in index_name:
        index_numbers.append(Lick_indices[index])
    return np.asarray(index_numbers)

#Return index reference number associated to a given Lick index
def I_number27(index_name):
    index_numbers=[]
    for index in index_name:
        index_numbers.append(Lick_indices27[index])
    return np.asarray(index_numbers)

#Return index reference number associated to a given Lick index
def I_number29(index_name):
    index_numbers=[]
    for index in index_name:
        index_numbers.append(Lick_indices29[index])
    return np.asarray(index_numbers)

#Return index reference number associated to a given Lick index
def I_number30(index_name):
    index_numbers=[]
    for index in index_name:
        index_numbers.append(Lick_indices30[index])
    return np.asarray(index_numbers)

#Return index reference number associated to a given Lick index
def I_numberall(index_name):
    index_numbers=[]
    for index in index_name:
        index_numbers.append(Lickemission_indices[index])
    return np.asarray(index_numbers)

#Return index reference number associated to a given Lick index
def I_numberopen(index_name):
    index_numbers=[]
    for index in index_name:
        index_numbers.append(Lickopen_indices[index])
    return np.asarray(index_numbers)

#Return index reference number associated to a given Lick index
def I_numberSDSS(index_name):
    index_numbers=[]
    for index in index_name:
        index_numbers.append(Lick_indicesSDSS[index])
    return np.asarray(index_numbers)

#Return index reference number associated to a given Lick index
def I_numberKnowles4(index_name):
    index_numbers=[]
    for index in index_name:
        index_numbers.append(Lick_Knowles4[index])
    return np.asarray(index_numbers)


#Return index reference number associated to a given Lick index
def I_typeopen(index_name):
    index_types=[]
    for index in index_name:
        index_types.append(Lickopen_type[index])
    return np.asarray(index_types)







#Returns the array of indices from a row of a Table
def array_de_table(ROW):
    return np.asarray([float(i) for i in ROW])



########################################################################
########################################################################
########################################################################
#Define the arrays of age, metallicity and alpha-enhancement
#parametros = np.load('./model/sMILES_velDisps/base/tabla/sMILES_base30_000_tabla.npy')
#
#t_array=np.unique(parametros[:, 0])
#Z_array=np.unique(parametros[:, 1])
#afe_array=np.unique(parametros[:, 2])
########################################################################
########################################################################
########################################################################






def I_model(t,Z,alpha,index_number, modelo, t_array, Z_array, afe_array):
    #Find the values of t, Z and alpha CLOSEST to the ones where the indices are already precomputed in the array of the model
    i=(np.abs(t_array-t)).argmin()
    j=(np.abs(Z_array-Z)).argmin()
    k=(np.abs(afe_array-alpha)).argmin()
    
    I=modelo[i,j,k,index_number]
    return I
    
def I_model_metal(Z,alpha,index_number, modelo, Z_array, afe_array):
    #Find the values of Z and alpha CLOSEST to the ones where the indices are already precomputed in the array of the model
    j=(np.abs(Z_array-Z)).argmin()
    k=(np.abs(afe_array-alpha)).argmin()
    
    I=modelo[j,k,index_number]
    return I


def get_data(measured_indices, index_numberlist):
    #index_numberlist is the ARRAY of index numbers to be used in the MCMC
    #data has shape (n,3), where n is the number of indices to be used in the MCMC: first column -> index values, second column -> index errors, third column -> index number
    data=measured_indices[index_numberlist,:]
    return data


#Define log-uniform distribution
def ln_uniform(value,a,b):
    if a<value<b:
        return 0
    else:
        return -np.inf

#Return the range of the uniform distribution of the prior of the parameter
def fpriors(variable, t_array, Z_array, afe_array):
    if variable=="t":
        return np.amin(t_array), np.amax(t_array)
    if variable=="Z":
        return np.amin(Z_array), np.amax(Z_array)
    if variable=="alpha":
        return np.amin(afe_array), np.amax(afe_array)

#Define the log-prior
def ln_prior(vals,var_nms, t_array, Z_array, afe_array):
    #vals=(t    index_number=math.floor(),Z,alpha) #var_nms=["t","Z","alpha"]
    l=0
    for value,variable in zip(vals,var_nms):
        l+=ln_uniform(value,fpriors(variable, t_array, Z_array, afe_array)[0],fpriors(variable, t_array, Z_array, afe_array)[1])
    return l

#Define the log-likelihood
def ln_like(vals,data, modelo, t_array, Z_array, afe_array):
    #vals=(t,Z,alpha)
    #data has shape (n,3), where n is the number of indices used
    n=len(data[:,0])
    I=data[:,0] #indices
    sigma=data[:,1] #errors 
    index_number=data[:,2]#number referencing each index
    
    index_number=index_number.astype(int) #CONVERT TO INTEGERS!
    
    Imodel=I_model(vals[0],vals[1],vals[2],index_number, modelo, t_array, Z_array, afe_array) #index according to model
    
    uno = -0.5*np.sum(((I-Imodel)/sigma)**2.0)
    dos = -1.0*np.sum(np.log(sigma))
    tres = -1.0*n*np.log((2.0*np.pi)**0.5)
    
    return uno + dos + tres
#return -0.5(np.sum(2.0*np.pi*sigma**2.0)+np.sum((I-Imodel)**2.0/sigma**2.0))
    

#Define the log-posterior
def ln_post(vals,var_nms,data, modelo, t_array, Z_array, afe_array):
    #vals=(t,Z,alpha) #var_nms=["t","Z","alpha"]
    
    #Check that the values are within the priors. Otherwise, return -Inf
    lpr=ln_prior(vals,var_nms, t_array, Z_array, afe_array)
    if not np.isfinite(lpr):
        return -np.inf
    #Compute the log-likelihood
    llh=ln_like(vals,data, modelo, t_array, Z_array, afe_array)
    
    return lpr+llh
    






#ONLY METAL FIT





#Return the range of the uniform distribution of the prior of the parameter
def fpriors_metal(variable, Z_array, afe_array):
    if variable=="Z":
        return np.amin(Z_array), np.amax(Z_array)
    if variable=="alpha":
        return np.amin(afe_array), np.amax(afe_array)

#Define the log-prior
def ln_prior_metal(vals,var_nms, Z_array, afe_array):
    l=0
    for value,variable in zip(vals,var_nms):
        l+=ln_uniform(value,fpriors_metal(variable, Z_array, afe_array)[0],fpriors_metal(variable, Z_array, afe_array)[1])
    return l

#Define the log-likelihood
def ln_like_metal(vals,data, modelo, Z_array, afe_array):
    n=len(data[:,0])
    I=data[:,0] #indices
    sigma=data[:,1] #errors 
    index_number=data[:,2]#number referencing each index
    
    index_number=index_number.astype(int) #CONVERT TO INTEGERS!
    
    Imodel=I_model_metal(vals[0],vals[1], index_number, modelo, Z_array, afe_array) #index according to model
    
    uno = -0.5*np.sum(((I-Imodel)/sigma)**2.0)
    dos = -1.0*np.sum(np.log(sigma))
    tres = -1.0*n*np.log((2.0*np.pi)**0.5)
    
    return uno + dos + tres
#return -0.5(np.sum(2.0*np.pi*sigma**2.0)+np.sum((I-Imodel)**2.0/sigma**2.0))
    

#Define the log-posterior
def ln_post_metal(vals,var_nms,data, modelo, Z_array, afe_array):
    #vals=(t,Z,alpha) #var_nms=["t","Z","alpha"]
    
    #Check that the values are within the priors. Otherwise, return -Inf
    lpr=ln_prior_metal(vals,var_nms, Z_array, afe_array)
    if not np.isfinite(lpr):
        return -np.inf
    #Compute the log-likelihood
    llh=ln_like_metal(vals,data, modelo, Z_array, afe_array)
    
    return lpr+llh


