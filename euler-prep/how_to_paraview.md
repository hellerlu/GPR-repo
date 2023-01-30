# Instruction on how to run paraview on Euler

## Introduction
Paraview can be run locally, which was found to be the most efficient and reliable. But if the models get larger, the RAM on local machines are not large enough to handle the task. Therefore, Euler was used. There are three ways to run Paraview on Euler:
- On the starting node with an interactive session and X11 forwarding:
    - Use the guide on [EulerDocs](https://scicomp.ethz.ch/wiki/ParaView) to start an interactive session 
- Only on Euler with pure python scripting. For this, run pvpython and a paraview script to get an output. More can be found on the [Paraview Wiki](https://www.paraview.org/Wiki/PvPython_and_PvBatch). 
- In a Client-Server-Mode using the pvpython implementations on Euler found on [Paraview Wiki](https://www.paraview.org/Wiki/PvPython_and_PvBatch), together with the instruction on [EulerDocs](https://scicomp.ethz.ch/wiki/ParaView_Client-Server) and the personal implementation described below.

## Personal implementation
I used Linux to interact with Euler in a Client-Server Mode, rather than Windows and MobaXterm. Ubuntu was installed and launched twice:
1. The first one connects to the server and runs the script `gprMax_pvScript.py` with pvpython. It then listens to connection on the predefined port in the script. 
2. The second one launches paraview locally and sends a request over the predefined port. 

As soon as the connection is established, the script continues and runs the model that is defined in the USER INPUT section. 

Obviously this only works with a valid VPN connection or in the ETH network. 

