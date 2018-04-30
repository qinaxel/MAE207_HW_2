import odrive.core
import time
import math

import numpy as np
import matplotlib.pyplot as plt

# For symbolic processing
import sympy
from sympy import symbols
from sympy import sin, cos, tan, asin, acos, atan, pi, diff
from sympy.utilities.lambdify import lambdify
from sympy import Matrix
from sympy import array
from sympy.solvers import solve

class Leg:
    """
    This is our first class in class :)

    We will define a leg class to interface with the leg and standardize 
    the kind of operations we want to perform

    """

    #### Variables outside the init function are constants of the class
    # leg geometry
    global l_base, l1, l2, encoder2angle, theta0_sym, theta1_sym, alpha0_sym, alpha1_sym
    
    l1 = 7  # NEED TO UPDATE units of cm
    l2 = 14  # NEED TO UPDATE units of cm
    l_base = 7.7  # NEED TO UPDATE units of cm
    theta0_sym, theta1_sym, alpha0_sym, alpha1_sym = symbols(
            'theta0_sym theta1_sym alpha0_sym alpha1_sym', real=True)
    # motor controller parameters
    encoder2angle = 2048 * 4

    ### Methods
    # Classes are initiated with a constructor that can take in initial parameters. At
    # a minimum it takes in a copy of itself (python... weird). The constructor
    # is one place we can define and initialize class variables

    def __init__(self, simulate = True):
        """
        This is the constructor for the leg class. Whenever you make a new leg
        this code will be called. We can optionally make a leg that will simulate
        the computations without needing to be connected to the ODrive
        """

        self.simulate = simulate

        # make the option to code without having the odrive connected
        if self.simulate == False:
            self.drv = self.connect_to_controller()
            self.m0 = self.drv.motor0  # easier handles to the motor commands
            self.m1 = self.drv.motor1
            
            m0_pos, m1_pos = self.get_joint_pos()
            self.joint_0_pos = m0_pos
            self.joint_1_pos = m1_pos

        else:
            self.drv = None
            self.joint_0_pos=pi/2
            self.joint_1_pos=pi/2

        # home angles
        self.joint_0_home = 0
        self.joint_1_home = 0
        # current positions
       

        # We will compute the jacobian and inverse just once in the class initialization.
        # This will be done symbolically so that we can use the inverse without having
        # to recompute it every time
        print('here2')
        self.J = self.compute_jacobian()
        # self.Jacobian_inv = self.Jacobian.inv

    def connect_to_controller(self):
        """
        Connects to the motor controller
        """
        drv = odrive.core.find_any(consider_usb=True, consider_serial=False)

        if drv is None:
            print('No controller found')
        else:
            print('Connected!')
        return drv
    
    
    ###
    ### Motion functions
    ###
     
    
    def get_joint_pos(self):
        """
        Get the current joint positions and store them in self.joint_0_pos and self.joint_1_pos in degrees.
        Also, return these positions using the return statement to terminate the function
        """
        # if simulating exit function
        if self.simulate == True:
            return (self.joint_0_pos, self.joint_1_pos)


        else: 
            self.joint_0_pos = self.m0.encoder.pll_pos/encoder2angle*(2*pi) +pi/2 
            self.joint_1_pos = self.m1.encoder.pll_pos/encoder2angle*(2*pi) +pi/2

        return (self.joint_0_pos, self.joint_1_pos)


    def set_home(self):
        """
        This function updates the home locations of the motors so that 
        all move commands we execute are relative to this location. 
        """
        # if simulating exit function
        if self.simulate == True:
            return

        else: 
            self.home = (self.m0.encoder.pll_pos/encoder2angle*(2*pi) +pi/2,self.m1.encoder.pll_pos/encoder2angle*(2*pi) +pi/2) 
        
        
        
           


    def set_joint_pos(self, theta0, theta1, vel0=0, vel1=0, curr0=0, curr1=0):
        """
        Set the joint positions in units of deg, and with respect to the joint homes.
        We have the option of passing the velocity and current feedforward terms.
        """
        # if simulating exit function
        if self.simulate == True:
            self.joint_0_pos = theta0
            self.joint_1_pos = theta1
        else: 
            self.get_joint_pos()
            self.m0.set_pos_setpoint((theta0-self.home[0])*encoder2angle/(2*pi)-pi/2 ,vel0,curr0) #encoder value = angle*(-4000)/pi determiend in first homework
            self.m1.set_pos_setpoint((theta1-self.home[1])*encoder2angle/(2*pi)-pi/2 ,vel1,curr1)


    def move_home(self):
        """
        Move the motors to the home position
        """
        # if simulating exit function
        if self.simulate == True:
            return

        else: # Your code here
            self.m0.set_pos_setpoint(self.home[0])
            self.m1.set_pos_setpoint(self.home[1])
        
        

    def set_foot_pos(self, x, y):
        """
        Move the foot to position x, y. This function will call the inverse kinematics 
        solver and then call set_joint_pos with the appropriate angles
        """
        # if simulating exit function
        if self.simulate == True:
           (theta_0, theta_1) = self.inverse_kinematics(x,y)
           self.set_joint_pos(theta_0,theta_1) 
           return (theta_0,theta_1) 
            
      
        else:
            (theta_0, theta_1) = self.inverse_kinematics(x,y)
            self.set_joint_pos(theta_0,theta_1) 
        
        
        

    def move_trajectory(self, tt, xx, yy):
        """
        Move the foot over a cyclic trajectory to positions xx, yy in time tt. 
        This will repeatedly call the set_foot_pos function to the new foot 
        location specified by the trajectory.
        """
        # if simulating exit function
        if self.simulate == True:
            theta0a=[]
            theta1a=[]
            alpha0a=[]
            alpha1a=[]
            for i in range(tt):
                (theta0,theta1)=self.inverse_kinematics(xx[i],yy[i])
                theta0a.append(theta0)
                theta1a.append(theta1)
                (alpha0,alpha1)=self.compute_internal_angles(theta0,theta1)
                alpha0a.append(alpha0)
                alpha1a.append(alpha1)
                
        
                
        else:
            for i in range(tt):
                self.set_foot_pos(xx[i],yy[i])
                print(theta_0,theta_1)
            
        return(theta0a,theta1a,alpha0a,alpha1a)
        

    ###
    ### Leg geometry functions
    ###
    def compute_internal_angles(self, theta_0, theta_1):
        """
        Return the internal angles of the robot leg 
        from the current motor angles
        """
        alpha_0, alpha_1, A, B, C = symbols('alpha_0 alpha_1 A B C', real=True)

        d = sympy.sqrt(l_base**2 + l1**2 - 2*l_base*l1*cos(theta_0))
        beta = -asin(l1/d*sin(theta_0))
        A = sympy.simplify(2*l1*l2*cos(theta_1)+2*d*l2*cos(beta))
        B = sympy.simplify(2*l1*l2*sin(theta_1)+2*d*l2*sin(beta))
        C = sympy.simplify(-l1**2-d**2-2*d*l1*cos(theta_1-beta))
        alpha_1 = sympy.simplify(atan(B/A) + acos(C/sympy.sqrt(A**2+B**2)))
        alpha_0 = acos((l1*cos(theta_1-beta) + l2*cos(alpha_1-beta)+d)/l2) + beta
        
        return (alpha_0, alpha_1)

    def compute_jacobian(self):
        """
        This function implements the symbolic solution to the Jacobian.
        """

        # initiate the symbolic variables
        #
        # Your code here that solves J as a matrix
        #
        (alpha0_sym,alpha1_sym) = self.compute_internal_angles(theta0_sym, theta1_sym)
        x = l_base/2 + l1*cos(theta0_sym) + l2*cos(alpha0_sym)
        y = l1*sin(theta1_sym) + l2*sin(alpha1_sym)
        
        J=Matrix([[sympy.diff(x,theta0_sym), sympy.diff(x,theta1_sym)],[sympy.diff(y,theta0_sym),sympy.diff(y,theta1_sym)]])
        return J

    def inverse_kinematics(self, x, y):
        """
        This function will compute the required theta_0 and theta_1 angles to position the 
        foot to the point x, y. We will use an iterative solver to determine the angles.
        """
        J = self.compute_jacobian()
        
        (theta_0,theta_1) = self.get_joint_pos() 
        #j_current = J.subs({theta0_sym: theta_current[0], theta1_sym: theta_current[1]})
        
        #J_inv=j_current.inv()
        
        
        x_target=Matrix([[x],[y]])
       
        #alpha=self.compute_internal_angles(theta_0,theta_1)
        #x1=l_base*0.5+l1*cos(theta_0)+l2*cos(alpha[0])
        #y1=l1*sin(theta_1)+l2*sin(alpha[1])
        #x_current=Matrix([[x1],[y1]])
        #x_error = Matrix([1e5, 1e5])
        # solution parameters
        # step size
        epsilon = 1e-2      # stop error
        
        for kk in range(1000):
            
            alpha=self.compute_internal_angles(theta_0,theta_1)
            x1=l_base*0.5+l1*cos(theta_1)+l2*cos(alpha[1])
            y1=l1*sin(theta_0)+l2*sin(alpha[0])
            x_current=Matrix([[x1],[y1]])
            x_error = (x_target - x_current)
            if x_error.norm() <  epsilon:
                break
            j_current=J.subs({theta0_sym: theta_0, 
                              theta1_sym: theta_1})
            J_inv=j_current.inv()
            d_theta = J_inv@x_error*1e-1 #how much to move theta
            theta_0 = theta_0 + d_theta[0] #update theta 
            theta_1 = theta_1 + d_theta[1]
            
        
        return (theta_0, theta_1)

    ###
    ### Visualization functions
    ###
    def draw_leg(self,l0, l1, l2, theta1, theta2, alpha1, alpha2, ax=False):
        """
        This function takes in the four angles of the leg and draws
        the configuration
        """

        link1, link2, width = l1, l2, l0
    
        if ax == False:
            ax = plt.gca()
            ax.cla()

        ax.plot(-width / 2, 0, 'ok')
        ax.plot(width / 2, 0, 'ok')
        
        ax.plot([-width / 2, 0], [0, 0], 'k')
        ax.plot([width / 2, 0], [0, 0], 'k')
        
        ax.plot(-width / 2 + np.array([0, link1 * cos(theta1)]), [0, link1 * sin(theta1)], 'k')
        ax.plot(width / 2 + np.array([0, link1 * cos(theta2)]), [0, link1 * sin(theta2)], 'k')
        
        ax.plot(-width / 2 + link1 * cos(theta1) + np.array([0, link2 * cos(alpha1)]), \
                link1 * sin(theta1) + np.array([0, link2 * sin(alpha1)]), 'k');
        ax.plot(width / 2 + link1 * cos(theta2) + np.array([0, link2 * cos(alpha2)]), \
                        np.array(link1 * sin(theta2) + np.array([0, link2 * sin(alpha2)])), 'k');
                        
        x_current = float(width / 2 + link1 * cos(theta2) + link2 * cos(alpha2))
        y_current = float(np.array(link1 * sin(theta2) + link2 * sin(alpha2)))
                        
        ax.plot(x_current, y_current, 'ro');
                        
        ax.axis([-(l1 + l2), l1 + l2, -l1, l1 + l2])
        ax.invert_yaxis()