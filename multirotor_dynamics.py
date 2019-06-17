'''
  Abstract Python class for multirotor dynamics

  You implementing class should define the methods u2(), u3(), u4()

  Based on:

   @inproceedings{DBLP:conf/icra/BouabdallahMS04,
     author    = {Samir Bouabdallah and Pierpaolo Murrieri and Roland Siegwart},
     title     = {Design and Control of an Indoor Micro Quadrotor},
     booktitle = {Proceedings of the 2004 {IEEE} International Conference on Robotics and 
                  Automation, {ICRA} 2004, April 26 - May 1, 2004, New Orleans, LA, {USA}},
      pages     = {4393--4398},
      year      = {2004},
      crossref  = {DBLP:conf/icra/2004},
      url       = {https:#doi.org/10.1109/ROBOT.2004.1302409},
      doi       = {10.1109/ROBOT.2004.1302409},
      timestamp = {Sun, 04 Jun 2017 01:00:00 +0200},
      biburl    = {https:#dblp.org/rec/bib/conf/icra/BouabdallahMS04},
      bibsource = {dblp computer science bibliography, https:#dblp.org}
    }
  
  Copyright (C) 2019 Simon D. Levy
 
  MIT License
'''

import numpy as np

class Dynamics(object):

    _STATE_X         = 0
    _STATE_X_DOT     = 1
    _STATE_Y         = 2
    _STATE_Y_DOT     = 3
    _STATE_Z         = 4
    _STATE_Z_DOT     = 5
    _STATE_PHI       = 6
    _STATE_PHI_DOT   = 7
    _STATE_THETA     = 8
    _STATE_THETA_DOT = 9
    _STATE_PSI       = 10
    _STATE_PSI_DOT   = 11

     # Might want to allow G to vary based on altitude
    _G  = 9.80665 

    
    def bodyZToInertial(bodyZ, rotation):
        '''
        Converts body frame to inertial frame using rotation angles
        Optimized for body X=Y=0
        '''
 
        phi   = rotation[0]
        theta = rotation[1]
        psi   = rotation[2]

        cph = np.cos(phi)
        sph = np.sin(phi)
        cth = np.cos(theta)
        sth = np.sin(theta)
        cps = np.cos(psi)
        sps = np.sin(psi)

        # This is the rightmost column of the body-to-inertial rotation matrix
        R = np.array([sph*sps + cph*cps*sth, cph*sps*sth - cps*sph, cph*cth])

        return bodyZ * R

    def inertialToBody(inertial, rotation):
        '''
        Converts inertiai frame to body frame using rotation angles

        See Section 5 of http://www.chrobotics.com/library/understanding-euler-angles
        '''
 
        phi   = rotation[0]
        theta = rotation[1]
        psi   = rotation[2]

        cph = np.cos(phi)
        sph = np.sin(phi)
        cth = np.cos(theta)
        sth = np.sin(theta)
        cps = np.cos(psi)
        sps = np.sin(psi)

        R = np.array([[cps*cth,                cth*sps,                   -sth], 
                      [cps*sph*sth - cph*sps,  cph*cps + sph*sps*sth,  cth*sph], 
                      [sph*sps + cph*cps*sth,  cph*sps*sth - cps*sph,  cph*cth]])

        return np.dot(R, inertial)

    def eulerToQuaternion(eulerAngles):
        '''
        Converts Euler angles to quaternion

        See https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
        '''
        # Convenient renaming
        phi = eulerAngles[0] / 2
        the = eulerAngles[1] / 2
        psi = eulerAngles[2] / 2

        # Pre-computation
        cph = np.cos(phi)
        cth = np.cos(the)
        cps = np.cos(psi)
        sph = np.sin(phi)
        sth = np.sin(the)
        sps = np.sin(psi)

        # Conversion
        return (
                cph * cth * cps + sph * sth * sps,
                cph * sth * sps - sph * cth * cps, 
                -cph * sth * cps - sph * cth * sps,
                cph * cth * sps - sph * sth * cps)

    def __init__(self, motorCount, b, d, m, l, Ix, Iy, Iz, Jr, maxrpm):
        '''
        Constructor accepts motor count, physical constants from article, and maximum RPMs
        '''

        self.motorCount = motorCount
        self.x = np.zeros(12)
        self.omegas = np.zeros(motorCount)

        self.b = b
        self.d = d
        self.m = m
        self.l = l
        self.Ix = Ix
        self.Iy = Iy
        self.Iz = Iz
        self.Jr = Jr

        self.maxrpm = maxrpm

        # Values computed in Equation 6
        self.U1 = 0     # total thrust
        self.U2 = 0     # roll thrust right
        self.U3 = 0     # pitch thrust forward
        self.U4 = 0     # yaw thrust clockwise
        self.Omega = 0  # torque clockwise

        # Radians per second for each motor
        self.omegas = np.zeros(motorCount)

        # Inertial-frame acceleration
        self.inertialAccel = np.zeros(3)

        # Flag for whether we're airborne
        self.airborne = False

        # Takeoff altitude, for detecting a crash
        self.zstart = 0

    def start(self, location, rotation, airborne=False):
        '''
        Initializes kinematic pose, with flag for whether we're airbone (helps with testing gravity).
        location = X,Y,Z
        rotation = phi, theta, psi
        airborne allows us to start on the ground (default) or in the air (e.g., gravity test)
        '''

        # Initialize state
        self.x[Dynamics._STATE_X]         = location[0]
        self.x[Dynamics._STATE_X_DOT]     = 0
        self.x[Dynamics._STATE_Y]         = location[1]
        self.x[Dynamics._STATE_Y_DOT]     = 0
        self.x[Dynamics._STATE_Z]         = location[2]
        self.x[Dynamics._STATE_Z_DOT]     = 0
        self.x[Dynamics._STATE_PHI]       = rotation[0]
        self.x[Dynamics._STATE_PHI_DOT]   = 0
        self.x[Dynamics._STATE_THETA]     = rotation[1]
        self.x[Dynamics._STATE_THETA_DOT] = 0
        self.x[Dynamics._STATE_PSI]       = rotation[2]
        self.x[Dynamics._STATE_PSI_DOT]   = 0

        # Initialize inertial frame acceleration in NED coordinates
        self.inertialAccel = Dynamics.bodyZToInertial(-Dynamics._G, rotation)

        # We can start on the ground (default) or in the air
        self.airborne = airborne

        # Remember our altitude at takeoff
        self.zstart = location[2]

    def update(self, dt):
        '''
        Updates state.
        dt time in seconds since previous update
        '''
        # Use the current Euler angles to rotate the orthogonal thrust vector into the inertial frame.
        # Negate to use NED.
        euler = np.array([self.x[6], self.x[8], self.x[10]])
        ned = Dynamics.bodyZToInertial(-self.U1/self.m, euler)

        # We're airborne once net downward acceleration goes below zero
        netz = ned[2] + Dynamics.g
        if not self.airborne:
            self.airborne = netz < 0

        # Once airborne, we can update dynamics
        if self.airborne:

            # Make some useful abbreviations
            phidot = self.x[Dynamics._STATE_PHI_DOT]
            thedot = self.x[Dynamics._STATE_THETA_DOT]
            psidot = self.x[Dynamics._STATE_PSI_DOT]

            dxdt = np.array([

                # Equation 12: compute temporal first derivative of state.
                self.x[Dynamics._STATE_X_DOT],  
                ned[0],
                self.x[Dynamics._STATE_Y_DOT],
                ned[1],
                self.x[Dynamics._STATE_Z_DOT],
                netz,
                phidot,
                psidot*thedot*(self.Iy-self.Iz)/self.Ix - self.Jr/self.Ix*thedot*self.Omega + self.l/self.Ix*self.U2,
                thedot,
                -(psidot*phidot*(self.Iz-self.Ix)/self.Iy + self.Jr/self.Iy*phidot*self.Omega + self.l/self.Iy*self.U3), 
                psidot,
                thedot*phidot*(self.Ix-self.Iy)/self.Iz   + self.l/self.Iz*self.U4
           ]) 

            # Compute state as first temporal integral of first temporal derivative
            self.x += dt * dxdt

            # Once airborne, inertial-frame acceleration is same as NED acceleration
            self.inertialAccel = np.copy(ned)

    def setMotors(self, motorvals):
        '''
        Uses motor values to implement Equation 6.
        motorvals in interval [0,1]
        '''
    
        # Convert the  motor values to radians per second
        self.omegas = motorvals * self.maxrpm * np.pi / 30

        # Compute overall torque from omegas before squaring
        self.Omega = self.u4(self.omegas)

        # Overall thrust is sum of squared omegas
        self.omegas = self.omegas ** 2
        self.U1 = self.b * self.omegas

        # Use the squared Omegas to implement the rest of Eqn. 6
        self.U2 = self.b * self.u2(self.omegas)
        self.U3 = self.b * self.u3(self.omegas)
        self.U4 = self.d * self.u4(self.omegas)

    def getState(self):
        '''
        Gets current state
        Returns (location, rotation, state, crashed), where:
            state = (angularVel, bodyAccel, inertialVel, quaternion)
            crashed = True if crashed, False otherwise
        '''
        # Get most values directly from state vector
        angularVel  = self.x[Dynamics._STATE_PHI_DOT:Dynamics._STATE_PHI_DOT+5:2]
        inertialVel = self.x[Dynamics._STATE_X_DOT:Dynamics._STATE_X_DOT+5:2]
        rotation    = self.x[Dynamics._STATE_PHI:Dynamics._STATE_PHI+5:2]
        location    = self.x[Dynamics._STATE_X:Dynamics._STATE_X+5:2]

        # Convert inertial acceleration and velocity to body frame
        bodyAccel = Dynamics.inertialToBody(self.inertialAccel, rotation)

        # Convert Euler angles to quaternion
        quaternion = Dynamics.eulerToQuaternion(rotation)

        # Make pose and state tuples
        state = angularVel, bodyAccel, inertialVel, quaternion

        # If we're airborne, we've crashed if we fall below ground level
        crashed =  (location[2] > self.zstart) if self.airborne else False

        return location, rotation, state, crashed
