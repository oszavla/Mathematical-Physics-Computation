import sympy as sp

# 1. Define all symbolic physical variables
m1, m2, l1, l2, g = sp.symbols('m1 m2 l1 l2 g', positive=True)
th1, th2, dth1, dth2 = sp.symbols('theta_1 theta_2 dtheta_1 dtheta_2', real=True)

# 2. Define the Cartesian to Angular kinematic mappings
# P is the inner product of the two rod vectors
P = l1 * l2 * sp.cos(th1 - th2) 
y1 = -l1 * sp.cos(th1)
v1_sq = (l1 * dth1)**2
v_rel_sq = (l2 * dth2)**2

# 3. Construct the Matrix Equation M * Lambda = R
M = sp.Matrix([
    [2 * l1**2, -2 * P],
    [-2 * P / m1, 2 * (m1 + m2) * l2**2 / (m1 * m2)]
])

R = sp.Matrix([
    [m1 * (g * y1 - v1_sq)],
    [-v_rel_sq]
])

# 4. Solve the linear system for lambda_1 and lambda_2
# M.solve(R) is computationally safer than M.inv() * R
Lambda = M.solve(R)

# 5. Convert Lagrange multipliers to rod tensions
# T1 = -2 * l1 * lambda_1
# T2 = -2 * l2 * lambda_2
T1 = sp.trigsimp(sp.simplify(-2 * l1 * Lambda[0]))
T2 = sp.trigsimp(sp.simplify(-2 * l2 * Lambda[1]))

# 6. Display the exact analytic forms
print("--- Tension in the Upper Rod (T1) ---")
sp.pprint(T1)

print("\n--- Tension in the Lower Rod (T2) ---")
sp.pprint(T2)