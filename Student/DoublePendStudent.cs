//============================================================================
using System;

public partial class DoublePendSim : Simulator
{
    // physical parameters
    double L1;   // Length of rod 1
    double L2;   // Length of rod 2
    double m1;   // mass of pendulum 1
    double m2;   // mass of pendulum 2

    //------------------------------------------------------------------------
    // RHSFuncDoublePend:  Evaluates the right sides of the differential
    //                     equations for the double pendulum
    //------------------------------------------------------------------------
    private void RHSFuncDoublePend(double[] xx, double t, double[] ff)
    {
        double theta1 = xx[0];
        double theta2 = xx[1];
        double omega1 = xx[2];
        double omega2 = xx[3];

        // Evaluate right sides of differential equations of motion
        //********************************************************************
        // Students enter equations of motion below
        //********************************************************************

        double cosTheta2 = Math.Cos(theta2);
        double sinTheta2 = Math.Sin(theta2);

        double A = (m1+m2)*L1*L1;
        double B = m2*L1*L2*cosTheta2;
        double C = B;
        double D = m2*L2*L2;
        double det = A*D - B*C;

        double R1 = m2*L1*L2*omega2*omega2*sinTheta2 - 
            (m1+m2)*g*L1*Math.Sin(theta1);
        double R2 = -m2*L1*L2*omega1*omega1*sinTheta2 - 
            m2*g*L2*Math.Sin(theta1+theta2);

        ff[0] = omega1;   // time derivative of state theta1
        ff[1] = omega2 - omega1;   // time derivative of state theta2
        ff[2] = (D*R1 - B*R2)/det;   
        ff[3] = (-C*R1 + A*R2)/det;   // time derivative of state omega2
    }

    //******************************************************************
    // Students enter energy calculations here. 
    //******************************************************************
    // Kinetic energy ----------
    public double KineticEnergy
    {
        get{
            double theta1 = x[0];
            double theta2 = x[1];
            double u1 = x[2];
            double u2 = x[3];

            //########## YOU NEED TO CALCULATE THIS ###########
            return 0.0; 
        }
    }

    // Potential energy
    public double PotentialEnergy
    {
         get{
            double theta1 = x[0];
            double theta2 = x[1];

            //########## YOU NEED TO CALCULATE THIS ###########
            return 0.0; 
        }
    }
}