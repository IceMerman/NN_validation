from scipy.io import loadmat

# Load .mat data into Python
#mat = loadmat('Demanda/BaseDatosNE_31_10_2022New_England.mat')
# load database
mat = loadmat('db.mat')
labels = mat['issecure']
VM = mat['VM']
VA = mat['VA']
# Reshape sload half of the data is active power, the other one is reactive power
SLOAD = mat['SLOAD']