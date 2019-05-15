## System Identification

### About

System Identification is done on a TREX 550 flybarless helicopter using the parameters used for the Yamaha RMAX helicopter as done in Bernard Mettler's paper.
There are 6 methods that requires testing namely

- *Genetic Algorithms*- Testing done
- *Invasive Weed Optimization*- Testing done
- *Simulated Annealing*
- *Biogeography-based Optimization*
- *Cultural Algorithm*
- *Imperialist Collective Optimization*

We aim to minimize the error between the graphs of the actual data and the flight data. 

The major changes have been done in the sphere(or cost) function of each of the methods.

---
### How to run IWO
- Open `iwo.m` and `sphere.m` (Changes have been made only in that)
- Open `best.mat` This workspace contains inr-the input matrix of 4 inputs recorded along the time steps
                                         outr- 10 outputs excluding rfb,c and d because we do not have a way to measure them
                                         t- time steps ( required for using lsim function)
					 population - pre defined initial population because random initializations result in NaN. 
					 *In `sphere.m` which is a cost function we input the array of values for the 40 parameters and create a state space model with that. Now we send the rc input values(inr) to the model and simulate the output. we compare the simulated output with the actual output by using correlation coefficient and least squares (it needs to be refined). this function returns a value that must be minimised.*

- Run `iwo.m`
- Once iwo is done with 1000 iterations it stops
- Go to BestSol structure in the workspace and the best array of the 40 parameters can be found under the variable position
- Copy this array vertically in the variable in the workspace called popul
- Now run PopulCheck.m to visualize the model
-----------------------------------------------------------------------------------
The other methods also follow a similar structure

---

Creator - [Navaneeth Krishnan](https://www.facebook.com/princeniken)

Maintainer - [Saumya Kumaar](https://dronefreak.bitbucket.io/)
