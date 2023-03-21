# TR-DRS
> The aim is to determine the distribution of time of flight (DTOF) for different combinations of scattering and absorption coefficients using White Monte Carlo (WMC) and to investigate the sensitivity of detection of grey matter tissue at different time intervals of light collection.

---
## Simulation steps
1. Place the segmented model and the position/direction of the probes in the `models` folder.
    ![image](https://user-images.githubusercontent.com/56038738/226594425-0830a037-bfac-497d-a3a7-0e0d7524bf9c.png)
2. Use `S1_make_the_sim_setting.m` to make the settings, including
    1. how many photons to simulate.
    2. how many SDSs to simulate.
    3. which mus combination to simulate.

    We only need to set the mus to simulate because the forward uses WMC, so the mua can be any combination.
    For example, the mus for each layer are set as follows:
    ```matlab
    layer_mus = {[50:25:350], [50:37.5:350], [10 19 28 37], [50:60:350]};     % mus for each layer, 1/cm
    ```
    Then there will be 13⨉9⨉4⨉6=2808 sets of mus combinations. And the program will automatically generate these combinations for you.
3. Use `thisPC_sim_wl_index.txt` to set the beginning and endding index for wavelength to simulate by this PC.
    ```
    % The beginning and endding index for wavelength to simulate by this PC
    1 2808
    ```
4. Use `S2_run_script.m` to run the simulation.
    You can set it to run many simulations one after the other.
    ![image](https://user-images.githubusercontent.com/56038738/226598088-12471065-8861-48fd-8d5b-dc410a4f5af2.png)
    If you have more than one GPU on your computer, you should set the `GPU_setting.txt` to determine which GPU is used and the load for each GPU.
    ```
    % Bool whether to use each GPU or not, 1*nGPU(amount of GPU) can be used by this machine
    1 1 1 0
    % workload for each GPU, 1*nGPU(amount of GPU) can be used by this machine
    105 109 100 25
    ```
5. The result folder for each subject will contain many [sim_ + index] folders
    ![image](https://user-images.githubusercontent.com/56038738/226600382-89e74761-f26b-4b99-89cb-591845672643.png)
    Each folder contains the WMC simulation result of the given MUS combination. It records the path length in each layer for each photon detected.
    There will also be some files containing the lookup table information, e.g. mus_table.txt containing the mus for each set of combinations.
    ![image](https://user-images.githubusercontent.com/56038738/226600841-b4c39bb7-4d64-407e-a4fb-d0b6402fd022.png)

---
## Calculate sensitivity
1. After simulation, use `S0_find_sim_prop.m` to select the layered optical properties (mua, mus, g, n) for given wavelengths.
    Here we use three wavelengths of 760, 806 and 850 nm.
2. Use `S1_make_the_sim_setting.m` to calculate the reflectance and pathlength from MC simulation, and split it into 25 timegates.
    ```matlab
    %% fun_MCX_sim_dist2axis.m
    cfg.tstart = 0;         % start of the simulation time window (in second)
    cfg.tend = 5e-9;        % end of the simulation time window (in second)
    cfg.tstep = 2e-10;      % time gate width (in second)

    %% S1_make_the_sim_setting.m
    num_gate = 25;          % number of time gates
    ```
3. Use `S2_GM_sensitivity_PL.m` to calculate grey matter sensitivity as the ratio of grey matter pathlength to total brain average pathlength.
    <img src="https://user-images.githubusercontent.com/56038738/226636556-3abc171c-9247-44d0-87dd-0673554f1500.png" width="500"/>
4. Use `S3_calculate_CV.m` to calculate the of Coefficient of Variation (CV) at different total photon counts.
    <img src="https://user-images.githubusercontent.com/56038738/226639292-e8ac1464-852b-45b1-9382-a1378c4fc50f.png" width="500"/>
5. Use `S4_GM_sensitivity_ref.m` to calculate the change in intensity of the reflectance of different time gates as the grey matter scattering coefficient changes.
    <img src="https://user-images.githubusercontent.com/56038738/226642411-61079335-3abf-4293-8158-94b468d5f4da.png" width="500"/>
    The sensitivity of each SDS is calculated by averaging the different variations and substituting them into the above simulation results above.
    <img src="https://user-images.githubusercontent.com/56038738/226642628-bca9cc74-8df2-460b-a82d-d7e5bb5bc6b4.png" width="500"/>
