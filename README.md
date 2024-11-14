# Turbulent Channel Flow
## Overview
Turbulent flow between two parallel plates is called channel flow.

## Usage
1.Make empty folder named 「output」 in same directory <br>
2.Make initial flow field by using 「initial.f90」<br>
3.Run the main program (grid counts must be same with 「initial.f90」). <br>
4.Visualize calculation result in 「output」 by paraview <br>

## Nondimensionalized Governing Equations
![ccc](https://github.com/user-attachments/assets/33ba14f5-5a8b-4cc6-8af6-57cb0de9b254)
<br>
boundary condition : x, z : periodic, y : wall (non-slip) <br>
initial condition : fully developed turbulent flow is independent of initial condition

<br>
<br>



## Result
![channelflow](https://github.com/user-attachments/assets/c2ce3327-a38a-4335-a494-10d3aca55244)
<br>
![aaa](https://github.com/user-attachments/assets/3b986e81-1d57-4706-836a-71773443c9ac)
contour of Q (Second invariant of the velocity gradient tensor) 
<br>
<br>
Voltex structure is found a little at near wall at initial time, but this is not fully developed turbulent voltex.
More trial with high Ret, more grid counts and time steps are needed.

## Comment
The figure is resulted with 
grid counts : 64,64,64 <br>
Ret : 300
