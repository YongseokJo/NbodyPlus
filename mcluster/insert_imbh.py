import numpy as np

########Physical Constants########
c = 2.99792458*10**10 #[cm s-1]
G = 6.67408*(10**(-8)) #[cm3 g-1 s-2]
pc = 3.0856776*10**18 #[cm]
Msun = 1.9891*10**33 #[g]

yr = 365*24*60*60
Myr = 1e6*yr


# Read Parameter and Compute N-body Unit
Input = []
file = open("./astro.PAR", "r")
while True:
    line = file.readline().split()
    Input.append(line)
    if not line:
        break
file.close()

N = int(Input[1][0])
RBAR = float(Input[2][7])
ZMBAR = float(Input[2][8])

Mtotal = N*ZMBAR # Msun
VBAR = np.sqrt(G*Mtotal*Msun/RBAR/pc)/(1e5)
TBAR = RBAR*pc/(VBAR*1e5)/(Myr)

# fort.10 파일로부터 입자들의 정보를 읽어온다.
Particle = np.loadtxt("./astro.NBODY", usecols = range(7))
# Particle을 질량순으로 정렬해준다.
Particle = Particle[Particle[:,0].argsort()[::-1]]


# 입자들의 정보에서 단위를 복원한다.
mass = Particle[:,0]*Mtotal
x = Particle[:,1]*RBAR
y = Particle[:,2]*RBAR
z = Particle[:,3]*RBAR
vx = Particle[:,4]*VBAR
vy = Particle[:,5]*VBAR
vz = Particle[:,6]*VBAR

kstar = np.full(N, int(14))
for i in range(N):
    if (mass[i] < 10.0):
        kstar[i] = int(1)

New_Particle = np.array([mass/Mtotal, x/RBAR, y/RBAR, z/RBAR, vx/VBAR, vy/VBAR, vz/VBAR, kstar]).T #kstar
np.savetxt("dat.10", New_Particle, fmt = "%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%.6e\t%d")

# Replace Input
Input[0][1] = str(1000000.0)
Input[1][0] = str(N)
Input[2][5] = str(1000.00)
Input[2][7] = str(RBAR)
Input[2][8] = str(ZMBAR)
Input[3][0] = str(0)
Input[4][3] = str(0)
Input[5][2] = str(2)
Input[7][5] = str(1)
Input[7][6] = str(1)

# Define the file name
file_name = "param.inp.old"

# Open the file in write mode
with open(file_name, 'w') as f:
    # Iterate over each sublist in the data list
    for sublist in Input:
        # Join the elements of the sublist into a string separated by a comma
        line = ' '.join(sublist)
        # Write the line to the file followed by a newline character
        f.write(line + '\n')
