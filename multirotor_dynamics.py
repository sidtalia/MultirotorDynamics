'''
  Python class for multirotor dynamics
 
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

class MultirotorDynamics(object):

    STATE_X         = 0
    STATE_X_DOT     = 1
    STATE_Y         = 2
    STATE_Y_DOT     = 3
    STATE_Z         = 4
    STATE_Z_DOT     = 5
    STATE_PHI       = 6
    STATE_PHI_DOT   = 7
    STATE_THETA     = 8
    STATE_THETA_DOT = 9
    STATE_PSI       = 10
    STATE_PSI_DOT   = 11

     # Might want to allow G to vary based on altitude
    G  = 9.80665 

    # bodyToInertial method optimized for body X=Y=0
    def bodyZToInertial(bodyZ, rotation):

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

    def __init__(self, motorCount, b, d, m, l, Ix, Iy, Iz, Jr, maxrpm):

        self.motorCount = motorCount
        self._x = np.zeros(12)
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
        
    def start(self, pose, airborne=False):
        '''
        Initializes kinematic pose, with flag for whether we're airbone (helps with testing gravity).
        pose location X,Y,Z rotation phi,theta,psi
        airborne allows us to start on the ground (default) or in the air (e.g., gravity test)
        '''

        location, rotation = pose

        # Initialize state
        self.x[MultirotorDynamics.STATE_X]         = location[0]
        self.x[MultirotorDynamics.STATE_X_DOT]     = 0
        self.x[MultirotorDynamics.STATE_Y]         = location[1]
        self.x[MultirotorDynamics.STATE_Y_DOT]     = 0
        self.x[MultirotorDynamics.STATE_Z]         = location[2]
        self.x[MultirotorDynamics.STATE_Z_DOT]     = 0
        self.x[MultirotorDynamics.STATE_PHI]       = rotation[0]
        self.x[MultirotorDynamics.STATE_PHI_DOT]   = 0
        self.x[MultirotorDynamics.STATE_THETA]     = rotation[1]
        self.x[MultirotorDynamics.STATE_THETA_DOT] = 0
        self.x[MultirotorDynamics.STATE_PSI]       = rotation[2]
        self.x[MultirotorDynamics.STATE_PSI_DOT]   = 0

        # Initialize inertial frame acceleration in NED coordinates
        self.inertialAccel = MultirotorDynamics.bodyZToInertial(-MultirotorDynamics.G, rotation)

        # We can start on the ground (default) or in the air
        self.airborne = airborne

        # Remember our altitude at takeoff
        self.zstart = pose.location[2]

    def update(self, dt):
        '''
        Updates state.
        dt time in seconds np.since previous update
        '''
        # Use the current Euler angles to rotate the orthogonal thrust vector into the inertial frame.
        # Negate to use NED.
        euler = np.array([self.x[6], self.x[8], self.x[10]])
        ned = MultirotorDynamics.bodyZToInertial(-self.U1/self.m, euler)

        # We're airborne once net downward acceleration goes below zero
        netz = ned[2] + MultirotorDynamics.g
        if not self.airborne:
            self.airborne = netz < 0

        # Once airborne, we can update dynamics
        if self.airborne:

            # Make some useful abbreviations
            phidot = self.x[MultirotorDynamics.STATE_PHI_DOT]
            thedot = self.x[MultirotorDynamics.STATE_THETA_DOT]
            psidot = self.x[MultirotorDynamics.STATE_PSI_DOT]

            dxdt = np.array([

                # Equation 12: compute temporal first derivative of state.
                self.x[MultirotorDynamics.STATE_X_DOT],  
                ned[0],
                self.x[MultirotorDynamics.STATE_Y_DOT],
                ned[1],
                self.x[MultirotorDynamics.STATE_Z_DOT],
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

'''
    protected:

        # roll right
        virtual double u2(double * o) = 0

        # pitch forward
        virtual double u3(double * o) = 0

        # yaw cw
        virtual double u4(double * o) = 0

        # motor direction for animation
        virtual int8_t motorDirection(uint8_t i) = 0

    public:

        /**
         * Exported state representations
         */

        # Kinematics
        typedef struct {

            double location[3]
            double rotation[3] 

        } pose_t

        typedef struct {

            double angularVel[3] 
            double bodyAccel[3] 
            double bodyVel[3] 
            double inertialVel[3] 
            double quaternion[4] 

            pose_t pose

        } state_t

        /**
         *  Destructor
         */
        virtual ~MultirotorDynamics(void)
        {
            delete _omegas
        }

       }

        /*
         *  Gets current state
         *
         *  @param state data structure for state
         *
         *  @return true if crashed, false otherwise
         */
        bool getState(state_t & state)
        {

            # Get most values directly from state vector
            for (uint8_t i=0 i<3 ++i) {
                uint8_t ii = 2 * i
                state.angularVel[i]    = _x[STATE_PHI_DOT+ii]
                state.inertialVel[i]   = _x[STATE_X_DOT+ii]
                state.pose.rotation[i] = _x[STATE_PHI+ii]
                state.pose.location[i] = _x[STATE_X+ii]
            }

            # Convert inertial acceleration and velocity to body frame
            inertialToBody(_inertialAccel, state.pose.rotation, state.bodyAccel)

            # Convert Euler angles to quaternion
            eulerToQuaternion(state.pose.rotation, state.quaternion)

            # If we're airborne, we've crashed if we fall below ground level
            return self.airborne ? (state.pose.location[2] > _zstart) : false
        }

        /**
         *  Supports debugging
         *
         * @return message string that can be displayed by calling program (e.g., FlightManager)
         */
        char * getMessage(void)
        {
            return _message
        }

        /**
         *  Frame-of-reference conversion routines.
         *
         *  See Section 5 of http:#www.chrobotics.com/library/understanding-euler-angles
         */

        static void bodyToInertial(double body[3], const double rotation[3], double inertial[3])
        {
            double phi   = rotation[0]
            double theta = rotation[1]
            double psi   = rotation[2]

            double cph = np.cos(phi)
            double sph = np.sin(phi)
            double cth = np.cos(theta)
            double sth = np.sin(theta)
            double cps = np.cos(psi)
            double sps = np.sin(psi)

            double R[3][3] = { {cps*cth,  cps*sph*sth - cph*sps,  sph*sps + cph*cps*sth}, 
                               {cth*sps,  cph*cps + sph*sps*sth,  cph*sps*sth - cps*sph}, 
                               {-sth,     cth*sph,                cph*cth}               }

            dot(R, body, inertial)
        }

        static void inertialToBody(double inertial[3], const double rotation[3], double body[3])
        {
            double phi   = rotation[0]
            double theta = rotation[1]
            double psi   = rotation[2]

            double cph = np.cos(phi)
            double sph = np.sin(phi)
            double cth = np.cos(theta)
            double sth = np.sin(theta)
            double cps = np.cos(psi)
            double sps = np.sin(psi)

            double R[3][3] = { {cps*cth,                cth*sps,                   -sth}, 
                               {cps*sph*sth - cph*sps,  cph*cps + sph*sps*sth,  cth*sph}, 
                               {sph*sps + cph*cps*sth,  cph*sps*sth - cps*sph,  cph*cth} }

            dot(R, inertial, body)
        }

        /**
         * Converts Euler angles to quaterion.
         *
         * @param eulerAngles input
         * @param quaternion output
         */

        static void eulerToQuaternion(const double eulerAngles[3], double quaternion[4])
        {
            # Convenient renaming
            double phi = eulerAngles[0] / 2
            double the = eulerAngles[1] / 2
            double psi = eulerAngles[2] / 2

            # Pre-computation
            double cph = np.cos(phi)
            double cth = np.cos(the)
            double cps = np.cos(psi)
            double sph = np.sin(phi)
            double sth = np.sin(the)
            double sps = np.sin(psi)

            # Conversion
            quaternion[0] =  cph * cth * cps + sph * sth * sps
            quaternion[1] =  cph * sth * sps - sph * cth * cps 
            quaternion[2] = -cph * sth * cps - sph * cth * sps
            quaternion[3] =  cph * cth * sps - sph * sth * cps
        }

        /**
         * Accessor method
         */
        uint8_t motorCount(void)
        {
            return _motorCount
        }

        /**
         * Factory method
         */
        static MultirotorDynamics * create(void)

} # class MultirotorDynamics
'''
