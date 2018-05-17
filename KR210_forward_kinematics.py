from sympy import symbols, cos, sin, pi, simplify
from sympy.matrices import Matrix
import numpy as np

#########################################################################################
# Create symbols, remember slices do not include the end value.
#########################################################################################
alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')  # q = theta

#########################################################################################
# KR210 DH parameters
#########################################################################################
# @formatter:off
s = {alpha0:     0, a0:      0, d1: 0.750,
     alpha1: -pi/2, a1:  0.350, d2:     0, q2: q2 - pi/2,
     alpha2:     0, a2:  1.250, d3:     0,
     alpha3: -pi/2, a3: -0.054, d4: 1.500,
     alpha4:  pi/2, a4:      0, d5:     0,
     alpha5: -pi/2, a5:      0, d6:     0,
     alpha6:     0, a6:      0, d7: 0.303, q7: 0}
# @formatter:on


#########################################################################################
# Homogeneous transforms
#########################################################################################
# @formatter:off
T0_1 = Matrix([[             cos(q1),            -sin(q1),            0,              a0],
               [ sin(q1)*cos(alpha0), cos(q1)*cos(alpha0), -sin(alpha0), -sin(alpha0)*d1],
               [ sin(q1)*sin(alpha0), cos(q1)*sin(alpha0),  cos(alpha0),  cos(alpha0)*d1],
               [                   0,                   0,            0,               1]])
T0_1 = T0_1.subs(s)

T1_2 = Matrix([[             cos(q2),            -sin(q2),            0,              a1],
               [ sin(q2)*cos(alpha1), cos(q2)*cos(alpha1), -sin(alpha1), -sin(alpha1)*d2],
               [ sin(q2)*sin(alpha1), cos(q2)*sin(alpha1),  cos(alpha1),  cos(alpha1)*d2],
               [                   0,                   0,            0,               1]])
T1_2 = T1_2.subs(s)

T2_3 = Matrix([[             cos(q3),            -sin(q3),            0,              a2],
               [ sin(q3)*cos(alpha2), cos(q3)*cos(alpha2), -sin(alpha2), -sin(alpha2)*d3],
               [ sin(q3)*sin(alpha2), cos(q3)*sin(alpha2),  cos(alpha2),  cos(alpha2)*d3],
               [                   0,                   0,            0,               1]])
T2_3 = T2_3.subs(s)

T3_4 = Matrix([[             cos(q4),            -sin(q4),            0,              a3],
               [ sin(q4)*cos(alpha3), cos(q4)*cos(alpha3), -sin(alpha3), -sin(alpha3)*d4],
               [ sin(q4)*sin(alpha3), cos(q4)*sin(alpha3),  cos(alpha3),  cos(alpha3)*d4],
               [                   0,                   0,            0,               1]])
T3_4 = T3_4.subs(s)

T4_5 = Matrix([[             cos(q5),            -sin(q5),            0,              a4],
               [ sin(q5)*cos(alpha4), cos(q5)*cos(alpha4), -sin(alpha4), -sin(alpha4)*d5],
               [ sin(q5)*sin(alpha4), cos(q5)*sin(alpha4),  cos(alpha4),  cos(alpha4)*d5],
               [                   0,                   0,            0,               1]])
T4_5 = T4_5.subs(s)

T5_6 = Matrix([[             cos(q6),            -sin(q6),            0,              a5],
               [ sin(q6)*cos(alpha5), cos(q6)*cos(alpha5), -sin(alpha5), -sin(alpha5)*d6],
               [ sin(q6)*sin(alpha5), cos(q6)*sin(alpha5),  cos(alpha5),  cos(alpha5)*d6],
               [                   0,                   0,            0,               1]])
T5_6 = T5_6.subs(s)

T6_7 = Matrix([[             cos(q7),            -sin(q7),            0,              a6],
               [ sin(q7)*cos(alpha6), cos(q7)*cos(alpha6), -sin(alpha6), -sin(alpha6)*d7],
               [ sin(q7)*sin(alpha6), cos(q7)*sin(alpha6),  cos(alpha6),  cos(alpha6)*d7],
               [                   0,                   0,            0,               1]])
T6_7 = T6_7.subs(s)
# @formatter:on

# Transform from base link to end effector (Intrinsic)
# T0_7 = simplify(T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_7)
T0_7 = simplify(T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_7)
print("T0_7 are:")
print(T0_7)

#########################################################################################
# Apply corrections due to gripper link coordinate frame has different definition
# in URDF file and DH parameters
#########################################################################################
# @formatter:off
R_z = Matrix([[    cos(np.pi), -sin(np.pi),             0, 0],
              [    sin(np.pi),  cos(np.pi),             0, 0],
              [             0,           0,             1, 0],
              [             0,           0,             0, 1]])

R_y = Matrix([[ cos(-np.pi/2),           0, sin(-np.pi/2), 0],
              [             0,           1,             0, 0],
              [-sin(-np.pi/2),           0, cos(-np.pi/2), 0],
              [             0,           0,             0, 1]])

R_correction = simplify(R_z * R_y)
# @formatter:on

# Apply correction term
T0_7_corrected = simplify(T0_7 * R_correction)
print("T0_7_corrected are:")
print(T0_7_corrected)

#########################################################################################
# Evaluate end effector coordinates
#########################################################################################
# Conversions between deg and rad
rtd = 180. / pi  # radians to degrees
dtr = pi / 180.  # degrees to radians

print("The coordinates of T0_7_corrected are:")
print(T0_7_corrected.evalf(subs={q1: 0,
                                 q2: 0,
                                 q3: 0,
                                 q4: 0,
                                 q5: 0,
                                 q6: 0,
                                 q7: 0}))

# When q1-7 are all 0, result should be:
# Matrix([
# [                  1.0,                     0, 6.12323399573677e-17, 2.153],
# [-7.49879891330929e-33,                   1.0, 1.22464679914735e-16,     0],
# [-6.12323399573677e-17, -1.22464679914735e-16,                  1.0, 1.946],
# [                    0,                     0,                    0,   1.0]]
# where (x,y,z) = (2.153, 0, 1.946) matches the reading under Displays-TF-Frames-gripper-link-position in RViz.
