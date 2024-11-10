WHAT THE REPOSITORY CONTAINS:
This is a suite of codes for solving the 1D Euler equations. For those unfamiliar, these equations govern the flow of a frictionless gas in one dimension. 
The key feature of the 1D Euler equations is they allow for shock wave formation owing to their inherent non-linear nature, making them a simple set of equations that
still captures some of the key physics that occurs during processes like supersonic flight, atmospheric re-entry, and stellar supernovae.
In order to study the behavior of shock waves in a variety of situations, scientists and engineers often turn to numerical solutions of the Euler equations.
However, properly simulating shock waves is a difficult task, particularly if one is interested in obtaining physically correct solutions with high orders of accuracy.
The critical difficulty lies in obtaining a proper numerical approximation for the flux difference occuring across the shock wave discontinuity, which if not done properly
can result in unphsyical solution oscillations (Gibbs phenomena), incorrect shock propagation speeds (convergence to the wrong weak solution), and numerical instability (unbounded growth of solution energy).
This repository contains several finite difference based numerical methods that accomplish the task of properly capturing shock waves to various degrees.
Several options for the spatial flux differencing scheme, time integration methods, and model problems for testing the numerical schemes are included.

SPATIAL FLUX SCHEMES
1) First order Roe scheme
2) WENO5 schemes as proposed by Prof. Chi Wang Shu and his collaborators
   
    i)   WENO5 Component wise reconstruction scheme
  
    ii)  WENO5 Characteristic wise reconstruction - Lax-Friedrichs flux splitting approach
    
    iii) WENO5 Characteristic wise reconstruction - Roe flux splitting approach
  
4) WENO6 schemes proposed by Hu and Adams
   
    i)   WENO6 Characteristic wise reconstruction - Lax-Friedrichs flux splitting approach

     ii)  WENO6 Characteristic wise reconstruction - Roe flux splitting approach

TIME INTEGRATION SCHEMES
1) First order explicit Euler method
2) Third-order TVD Runge-Kutta method ("SSP-RK3") of Shu and Osher

1D TEST PROBLEMS
1) Sod's shock tube problem
2) Shu-Osher entropy wave/shock interaction
3) A steady shock wave
4) The Lax problem
5) The 123 problem

SOLUTION VISUALIZATION AND COMPARISION
1) Built in plotting routines for 1D profiles
2) Ability to view x-t plane to visualze characteristic wave fronts
3) Exact solutions (computed on fine meshes w/ high order schemes) for the following problems
   
     i)  Sod problem
   
      ii) Shu-Osher problem


WHAT THE REPOSITORY IS FOR:
This repository started out as a project for a class I took on numerical methods for hyperbolic conservation laws, in order to help me better understand concepts
like the Roe scheme, WENO schemes, and the more general principles of turning the theory of hyperbolic PDE solution techniques into actual working code.
As such, the code in this repository is designed to be educational, and is written with copious comments that describe the methods being used and the papers and textbooks 
that are useful in understanding them. This might be especially useful for the new graduate student, scientist, or engineer who has not been directly exposed to the concepts
of hyperbolic PDE theory, applied math, or computational gas dynamics before and is keen on learning how to implement them.

The code is also written in a highly modular fashion, with arguably more subroutines than necessary. This was done intentionally in the hopes of partitioning the code into digestable
chunks that corresponde to a single step in the solution process (i.e, initial condition creation, flux reconstruction, time marching, plotting, etc.) and also to make the overall code 
easier to interact with.

I also want to stress that since the routines were written with education and user-friendliness in mind, they are not optimized for performance. For 1D problems, this really does not matter
but if one wants to extend these concepts to effective 2D or 3D solvers, code optimization would be required.  




HOW TO INTERACT WITH THE CODE:

The entirety of the code is controlled by the routine 'DRIVERSCRIPT.m'.
Within this routine are the options to change the flux scheme, the time integration scheme, the test problem, and the visualization approaches.
Comments are contained within this script that actually detail how to accomplish each of these tasks.
Furthermore, each subroutine that the code calls contains comments that describe its purpose and document the mathematical steps occurring in
each numerical algorithm. References to usefule papers and textbooks are included where appropriate for the user to find more detail when desired.
While the published code is written in Matlab, hopefully conversion to a different language can be accomplished in a relatively straightforward manner
through the use of some online GPT.


A HANDFUL OF FIGURES:

Herein I've included several figures that demonstrate the validity of the code.


Verification of 6th order accuract for WENO6 scheme - done by comparing a smooth numerical solution to Burgers' equation to an exact solution.
Exact solution obtained through Newtons method prior to shock formation. Order of accuracy is computed using the L1-norm of the numerical solution.

![burgers_weno6](https://github.com/sulli72/1D_EULER/assets/37673021/b26c0c2b-03f8-4acf-a544-491eaafc9853)

Copmuted solution to the Sod shock tube problem using WENO5 - Roe formulation. Shown are the temporally developing 1D variable profiles within the domain. 
The second figure shows the x-t plane for density, highlighting the expansion fan, contact wave, and shock wave that all develop
from the discontinuity in the initial condition.

![sod_soln](https://github.com/user-attachments/assets/7c1650e0-e177-47f6-af5d-892a198bb6be)



![density_chars_sod](https://github.com/sulli72/1D_EULER/assets/37673021/813b508a-aba7-4398-9952-d6d7051e1cac)

Finally, here is a comparison between the WENO5 and WENO6 methods (both using the Roe approach to wave speed determination). Since the WENO6 scheme
retains a symmetric stencil in smooth regions, it captures the high wavenumber features of the entropy disturbance significantly better than the
dissipative WENO5 scheme using the same mesh. Note that the 'exact' solution was computed on a 4000 point mesh using the WENO5 scheme. 
Also included is the time-varying solution using the WENO6 scheme.



i)  WENO5
![WENO5-SHU-V2](https://github.com/sulli72/1D_EULER/assets/37673021/22d4dff9-9faf-4de9-befd-27c935b13542)

ii) WENO6
![WENO6-SHU-V2](https://github.com/sulli72/1D_EULER/assets/37673021/87f1363f-1976-437c-8b45-759d426f0412)

iii) WENO6 SOLUTION 

![shu_soln](https://github.com/user-attachments/assets/56e27206-9a91-42ae-9239-02fd715457c1)

   
