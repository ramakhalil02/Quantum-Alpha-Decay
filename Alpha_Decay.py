import numpy as np
import matplotlib.pyplot as plt


#Constants 
hbar_c = 197.327        # MeV·fm
alpha = 1 / 137.036
amu_to_MeV = 931.494 
c = 299792458       #m / s

# Constants for 212Po 
Z_alpha = 2
Z_daughter = 82
A_alpha = 4
A_daughter = 208
Q = 8.954  # MeV
R_atomic = 9.01
R_max = 26.374
N = 4000

# Constants for 238U
# Z_alpha = 2
# Z_daughter = 90        
# A_alpha = 4
# A_daughter = 234       
# Q = 4.2749             # MeV
# R_atomic = 9.45        
# R_max = 60.63           


def reduced_mass(A1, A2):  # MeV/c^2
    return (A1 * A2) / (A1 + A2) * amu_to_MeV  

def coulomb_potential(r, Z1, Z2): #MeV
    return alpha * Z1 * Z2 * hbar_c / r

def compute_k(mu, Q, V):
    V0 = -25.0   #Depth of nuclear well
    V = np.insert(V, 0, V0)
    k = np.zeros_like(V, dtype=complex)
    for i in range(len(V)):
        if Q >= V[i]:
            k[i] = (1j * np.sqrt(2 * mu * (Q - V[i])) / hbar_c ) 
        else:
            k[i] = (np.sqrt(2 * mu * (V[i] - Q)) / hbar_c) 
    return k

def build_matrix(N, x, k):
    size = 2 * N + 3  
    M = np.zeros((size, size), dtype=complex)
    b = np.zeros(size, dtype=complex)

    # Incoming wave A = 1 
    b[0] = 1
    M[0, 0] = 1  # A

    #Outermost left boundary 
    
    M[1, 0] = np.exp(k[0] * x[0])     
    M[1, 1] = np.exp(-k[0] * x[0])   
    M[1, 2] = -np.exp(k[1] * x[0])    
    M[1, 3] = -np.exp(-k[1] * x[0])  

    M[2, 0] = k[0] * np.exp(k[0] * x[0])
    M[2, 1] = -k[0] * np.exp(-k[0] * x[0]) 
    M[2, 2] = -k[1] * np.exp(k[1] * x[0])
    M[2, 3] = k[1] * np.exp(-k[1] * x[0])

    # Internal boundaries 
    for i in range(1, N):
        xi = x[i]
        idx1 = 2 * i + 1
        idx2 = 2 * i + 2
        ai = 2 * i
        bi = 2 * i + 1
        aip1 = 2 * i + 2
        bip1 = 2 * i + 3

        e_l = np.exp(k[i] * xi)
        e_lm = np.exp(-k[i] * xi)
        e_r = np.exp(k[i+1] * xi)
        e_rm = np.exp(-k[i+1] * xi)

        M[idx1, ai] = e_l
        M[idx1, bi] = e_lm
        M[idx1, aip1] = -e_r
        M[idx1, bip1] = -e_rm

        M[idx2, ai] = k[i] * e_l
        M[idx2, bi] = -k[i] * e_lm
        M[idx2, aip1] = -k[i+1] * e_r
        M[idx2, bip1] = k[i+1] * e_rm

    # Last boundary 
    xf = x[-1]

    M[-2, -3] = np.exp(k[-2] * xf)    
    M[-2, -2] = np.exp(-k[-2] * xf)  
    M[-2, -1] = -np.exp(k[-1] * xf)  

    M[-1, -3] = k[-2] * np.exp(k[-2] * xf)
    M[-1, -2] = -k[-2] * np.exp(-k[-2] * xf)
    M[-1, -1] = -k[-1] * np.exp(k[-1] * xf)

    return M, b


def get_v_f(m, E):  #fm, s^-1
    v_by_c_2 = ( 2 * E ) / m
    v = np.sqrt(v_by_c_2) * c
    v_fm = v  * (10**15)
    f = v_fm / (2*R_atomic)
    return v_fm , f

def compute_lifetime(T, f):
    tau = 1 / (T * f)
    t_half = tau * np.log(2)
    t_half_ns = np.abs( t_half * 10**(9) )
    return tau, t_half_ns

def compute_wavefunction(solution, x, k, points_per_region=100):
    N = (len(solution) - 3) // 2
    A, B = solution[0], solution[1]
    F = solution[-1]

    r_full = []
    psi_full = []

    r_I = np.linspace(0, x[0], points_per_region)
    k0 = k[0]
    psi_I = A * np.exp(k0 * r_I) + B * np.exp(-k0 * r_I)
    r_full.extend(r_I)
    psi_full.extend(psi_I)

    for i in range(N):
        r_seg = np.linspace(x[i], x[i+1], points_per_region)
        a_i = solution[2 + 2 * i]
        b_i = solution[3 + 2 * i]
        psi_seg = a_i * np.exp(k[i+1] * r_seg) + b_i * np.exp(-k[i+1] * r_seg)
        r_full.extend(r_seg)
        psi_full.extend(psi_seg)

    r_III = np.linspace(x[-1], x[-1] + 10, points_per_region)
    kF = k[-1]
    psi_III = F * np.exp(kF * r_III)
    r_full.extend(r_III)
    psi_full.extend(psi_III)

    return np.array(r_full), np.array(psi_full)

def run_simulation(Q_value, N, R, R_max):
    mu = reduced_mass(A_alpha, A_daughter)
    x = np.linspace(R, R_max, N+1)
    dx = x[1] - x[0]
    V = coulomb_potential(x[:-1] + dx/2, Z_alpha, Z_daughter)
    k = compute_k(mu, Q, V)
    M, b = build_matrix(N, x, k)
    solution = np.linalg.solve(M, b)

    A = 1
    F = solution[-1]
    k_I = k[0]/(1j)
    k_III = k[-1]/(1j)
    T = abs(F / A)**2 * (k_III / k_I)

    v ,f = get_v_f((A_alpha*amu_to_MeV), Q)
    tau, t_half = compute_lifetime(T, f)
    
    r_vals, psi_vals = compute_wavefunction(solution, x, k)
    prob_density = np.abs(psi_vals)**2
    
    plt.figure(figsize=(10, 5))
    plt.plot(r_vals, np.log10(prob_density), label='Log Probability Density')
    plt.axvline(x=R, color='gray', linestyle='--', label='Boundary: Nucleus-Barrier')
    plt.axvline(x=R_max, color='gray', linestyle='--', label='Boundary: Barrier-Free region')
    plt.text(R - 4,-3, 'Region I', fontsize=10)
    plt.text((R + R_max)/2 - 3, -1, 'Region II', fontsize=10)
    plt.text(R_max + 1, -11, 'Region III', fontsize=10)
    
    plt.xlabel('Distance from nuclear center r [fm]')
    plt.ylabel(r'$\log_{10}(|\psi(r)|^2)$')
    plt.title('Log Probability Density of Alpha Particle')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

    return t_half


half_life_ns = run_simulation(Q, N, R_atomic, R_max)





      
      
