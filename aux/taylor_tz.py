# Vendored from the author's personal library (MYLIBS/TAYLOR_emcee_tz.py),
# verbatim -- the log-posterior for the cosmographic (Taylor-expansion)
# t(z) fit, narrow-priors variant. See aux/taylor_tz_wide.py for the wide
# variant (desicc/config.py's COSMOGRAPHIC_PRIORS_CHOICES, --priors in
# pipeline/30_cosmographic.py).

import numpy as np
import emcee
from scipy.integrate import simpson
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################


#Define log-uniform distribution
def ln_uniform(value,a,b):
    if a<value<b:
        return 0
    else:
        return -np.inf


#Return the range of the uniform distribution of the prior of the parameter
def fpriors(variable):
    if variable.startswith("age"):
        return 0.0, 20.0
    if variable=="Hz0":
        return 20.0, 220.0
    if variable=="qz0":
        return -2.0, 1.0
    if variable=="jz0":
        return -1.0, 3.0


#Define the log-prior
def ln_prior(vals,var_nms, orden):
    l = 0
    for value, variable in zip(vals, var_nms):
        l += ln_uniform(value, fpriors(variable)[0], fpriors(variable)[1])

    return l


# Nueva función lat que acepta un diccionario de edades presentes
def lat(orden, grupo, z, z0, ages_dict, Hz0, qz0, jz0):
    # 'grupo' es un array de strings como ['160180', '200225', ...]
    # Buscamos la edad correspondiente en el diccionario dinámico
    age0 = np.array([ages_dict[g] for g in grupo])

    constante = 978.5641806189751
    CtH = constante/Hz0
    y = (z-z0)/(1.0 + z0)

    C1 = 1.0 if orden >= 1 else 0
    C2 = -1.0*(2.0 + qz0)/2.0 if orden >= 2 else 0
    C3 = (1.0 + qz0 + 0.5*qz0**2.0 - 1.0/6.0 * jz0**2.0) if orden >= 3 else 0

    return age0 - CtH*(C1*y + C2*y**2.0 + C3*y**3.0)

def ln_like(orden, grupo, tobs, dtobs, z, z0, vals, var_nms):
    # Separamos parámetros de cosmología de los de edad
    # Asumimos que Hz0, qz0, jz0 son los últimos 3
    n_ages = len(vals) - 3
    ages_dict = dict(zip(var_nms[:n_ages], vals[:n_ages]))
    Hz0, qz0, jz0 = vals[-3:]

    # Limpiamos los nombres de las llaves del diccionario para que coincidan con 'grupo'
    # 'age160180' -> '160180'
    clean_ages_dict = {k.replace('age', ''): v for k, v in ages_dict.items()}

    tmodel = lat(orden, grupo, z, z0, clean_ages_dict, Hz0, qz0, jz0)

    chi2 = -0.5 * np.sum(((tmodel - tobs) / dtobs)**2.0)
    log_norm = -1.0 * len(tobs) * np.log(np.sqrt(2.0 * np.pi)) - np.sum(np.log(dtobs))
    return chi2 + log_norm



#Define the log-posterior
def ln_post(vals, var_nms, tobs, dtobs, z, z0, grupo, orden=3):
    lpr = ln_prior(vals, var_nms, orden)
    if not np.isfinite(lpr):
        return -np.inf

    # EL CAMBIO CLAVE: Pasa 'vals' y 'var_nms' como objetos, NO expandidos con *
    llh = ln_like(orden, grupo, tobs, dtobs, z, z0, vals, var_nms)

    return lpr + llh
