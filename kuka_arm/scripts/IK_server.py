#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np


class ForwardKinematicsSolver:
    """
    Solvers forward kinematics equations upon starting the server.
    """

    def __init__(self):
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
        # 1: homogeneous transforms between bas_link and gripper_link (end effector)
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

        # Transform from base_link to other joints
        T0_3 = simplify(T0_1 * T1_2 * T2_3)

        # Transform from base link to gripper_link (end effector) - intrinsic
        # T0_7 = simplify(T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_7)
        T0_7 = simplify(T0_1 * T1_2 * T2_3 * T3_4 * T4_5 * T5_6 * T6_7)

        #########################################################################################
        # 2: roll, pitch and yaw rotation between base_link and gripper_link
        #########################################################################################
        roll, pitch, yaw = symbols("roll pitch yaw")

        # @formatter:off
        # Rotation about X
        R_X = Matrix([[          1,           0,          0],
                      [          0,   cos(roll), -sin(roll)],
                      [          0,   sin(roll), cos(roll)]])
        # Rotation about Y
        R_Y = Matrix([[ cos(pitch),           0, sin(pitch)],
                      [          0,           1,          0],
                      [-sin(pitch),           0, cos(pitch)]])
        # Rotation about Z
        R_Z = Matrix([[   cos(yaw),   -sin(yaw),          0],
                      [   sin(yaw),    cos(yaw),          0],
                      [          0,           0,          1]])
        # @formatter:on

        R0_7 = R_Z * R_Y * R_X

        #########################################################################################
        # Apply corrections due to gripper link coordinate frame has different definition
        # in URDF file and DH parameters
        #########################################################################################
        # @formatter:off
        T_z = Matrix([[    cos(np.pi), -sin(np.pi),             0, 0],
                      [    sin(np.pi),  cos(np.pi),             0, 0],
                      [             0,           0,             1, 0],
                      [             0,           0,             0, 1]])

        T_y = Matrix([[ cos(-np.pi/2),           0, sin(-np.pi/2), 0],
                      [             0,           1,             0, 0],
                      [-sin(-np.pi/2),           0, cos(-np.pi/2), 0],
                      [             0,           0,             0, 1]])

        T_correction = simplify(T_z * T_y)
        R_correction = simplify(T_z[0:3, 0:3] * T_y[0:3, 0:3])
        # @formatter:on

        # Apply correction term
        T0_7_corrected = simplify(T0_7 * T_correction)
        R0_7_corrected = simplify(R0_7 * R_correction)

        # Put together member variables
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3
        self.q4 = q4
        self.q5 = q5
        self.q6 = q6
        self.q7 = q7
        self.T0_3 = T0_3
        self.T0_7 = T0_7
        self.R0_7 = R0_7
        self.T_correction = T_correction
        self.R_correction = R_correction
        self.T0_7_corrected = T0_7_corrected
        self.R0_7_corrected = R0_7_corrected


SIDE_A = 1.501
SIDE_C = 1.25

print("init FK solver...")
_fks = ForwardKinematicsSolver()
print("FK solver init done")


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        # Initialize service response
        joint_trajectory_list = []
        for pose_i in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            # roll, pitch, yaw = end-effector orientation
            px = req.poses[pose_i].position.x
            py = req.poses[pose_i].position.y
            pz = req.poses[pose_i].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[pose_i].orientation.x, req.poses[pose_i].orientation.y,
                 req.poses[pose_i].orientation.z, req.poses[pose_i].orientation.w])

            # Rotation matrix from base frame to end effector
            R0_7_corrected = _fks.R0_7_corrected.subs({"roll": roll,
                                                       "pitch": pitch,
                                                       "yaw": yaw})

            # Knowing r_ee_0 - the final coordinate of the end effector.
            r_ee_0 = Matrix([[px],
                             [py],
                             [pz]])
            # print("r_ee_0 are:")
            # print(r_ee_0)

            # Vector from wrist center to gripper
            # @formatter:off
            r_ee_wc = 0.303 * R0_7_corrected*Matrix([[0],
                                                     [0],
                                                     [1]])
            # @formatter:on

            # Get wrist center vector
            # From vector composition equation where r_ee_0 = r_wc_0 + r_ee_wc, find r_wc_0 = r_ee_0 - r_ee_wc.
            # Here r_wc_0 is expressed in terms of px, py, pz
            # All expressed in base frame
            r_wc_0 = r_ee_0 - r_ee_wc
            # print("r_wc_0 are:")
            # print(r_wc_0)

            # Get joint angles using trigonometric
            # Theta angles naming per course material Project: Robotic Arm - 15
            theta1 = atan2(r_wc_0[1], r_wc_0[0])
            side_b = sqrt(pow(sqrt(r_wc_0[0] * r_wc_0[0] + r_wc_0[1] * r_wc_0[1]) - 0.35, 2) +
                          pow(r_wc_0[2] - 0.75, 2))

            angle_a = acos((side_b * side_b + SIDE_C * SIDE_C - SIDE_A * SIDE_A) / (2 * side_b * SIDE_C))
            angle_b = acos((SIDE_A * SIDE_A + SIDE_C * SIDE_C - side_b * side_b) / (2 * SIDE_A * SIDE_C))

            theta2 = pi / 2. - angle_a - atan2(r_wc_0[2] - 0.75,
                                               sqrt(r_wc_0[0] * r_wc_0[0] + r_wc_0[1] * r_wc_0[1]) - 0.35)

            theta3 = pi / 2. - (angle_b + 0.036)

            R0_3 = _fks.T0_3[0:3, 0:3]
            R0_3 = R0_3.evalf(subs={_fks.q1: theta1,
                                    _fks.q2: theta2,
                                    _fks.q3: theta3})
            R3_6 = R0_3.inv("LU") * R0_7_corrected

            theta4 = atan2(R3_6[2, 2], -R3_6[0, 2])
            theta5 = atan2(sqrt(R3_6[0, 2] * R3_6[0, 2] + R3_6[2, 2] * R3_6[2, 2]), R3_6[1, 2])
            theta6 = atan2(-R3_6[1, 1], R3_6[1, 0])

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()


if __name__ == "__main__":
    IK_server()
    efasef = 5
