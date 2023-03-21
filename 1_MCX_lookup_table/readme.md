# MCX lookup table
Do the MCX simulation to generate the lookup table.

---

## Prepare
* The models will be used to simulate (including the position and direction of the probes)  
* The mus points for each layer you want to simulate.  

---

## Simulation steps

1. Put the segmented model and position/direction of the probes in the `models` folder.  
    ![](https://i.imgur.com/BUZHzL1.png)  

2. Use the `S1_make_the_sim_setting.m` to make the setting, including:  
    1. how many SDS are there to simulate.
    2. how many photon in each simultion.
    3. which mus combination to simulation.

    We only need to set the mus to simulate, because the forward is using WMC, so the mua can be any combination.  
    For example, the mus for each layer are setting as below:  
    ![](https://i.imgur.com/OdeYvgx.png)  
    Then there will be 13X9X4X6=2808 sets of mus combinations. And the program will auto generate these combinations for you.  

3. Use `S2_run_script.m` to run the simulation.  
    You can set it to run many simulation one-by-one.  
    ![](https://i.imgur.com/Tdp4faW.png)  
    If you hany more than one GPU on your computer, you should set the `GPU_setting.txt` to determine which GPU is used and the load for each GPU.  
    ![](https://i.imgur.com/YEHCKdz.png)  

4. The result folder for each subject will containe many [sim_ + index] folders  
    ![](https://i.imgur.com/VeaoCSL.png)  
    In each folder is the WMC simulation result of the given mus combination. It records the pathlength in each layer for each detected photon.  
    There will also be some files containing the information of the lookup table, e.g., the `mus_table.txt` containing the mus for each set of combination.  
    ![](https://i.imgur.com/IY0smvV.png)  

5. After the simulation, use `S3_exam_the_simed_table.m` to check if there is any error in the simulation result.  

6. If you want to make additional mus combination, use `S4_make_additional_sim_setting.m` to do that, it will compare the new mus sets to the old one, and put the additional mus sets in the end of the mus list (so the original simulated result should not be change.)