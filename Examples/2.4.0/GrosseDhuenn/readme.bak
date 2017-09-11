# *Grosse Dhuenn* Reservoir

## Study Site
*Grosse Dhuenn* is a [multi-purpose reservoir](http://www.sciencedirect.com/science/article/pii/S0273122398000432) (51.07 ºN, 7.22 ºE) situated in a densely populated area of Cologne, Western Germany.
The reservoir has a surface area of 3.7501 km<sup>2</sup>, and is about 48 m deep.

## Model Setup
The model is set-up to simulate the hydrological domain of *Grosse Dhuenn* reservoir for 2 years, from **1996-01-01** to **1998-01-01** (`&time`), with water balance and heat fluxes [hypothetically calculated](http://aed.see.uwa.edu.au/research/models/GLM/downloads/AED_GLM_v2_0b0_20141025.pdf) based on the reservoir configuration and input data. Water quality functionality (AED2) is disabled for the simulation (`&wq_setup`).


#### Reservoir Configuration
The flexible model structure allow users to configure the switches of the individual model components, enabling the customization of the physical model, meteorological and flow conceptualisations.


**Physical Model Configuration**

| Configuration Flags | Description | Activation |
| ---------------- |:----------|:-----------:|
| `&morphometry` | Reservoir location (*latitude and longitude*) and bathymetry (*depth and area*) setup | - |
| `&wq_setup` | Water Quality setup | *Disabled* |
| `&init_profiles` | Initial reservoir conditions with depths (eg. *temperature and salinity* ) | - |
| `&output` | Specification of output file details | - |

**Meteorological Configuration**

| Configuration Flags | Description | Activation|
| ---------------- |:----------|:-----------:|
| `met_sw = .true.` | Surface meteorological forcing | *Enabled* |
| `subdaily = .true.` | Subdaily meteorological input setup switch (in *hourly timesteps*; `dt = 3600`) | *Enabled* |
| `&snowice` | Snow-ice lake dynamics setup | *Enabled* |
| `&bird_model` | Incorporation of *Bird Solar Radiation* model which theoretically computes solar radiation based on available data| *Enabled* |
-------How about rad_mode? compute longwave from cloudcover by default?

**Flow Configuration**

| Configuration Flags | Description | Activation |
| ---------------- |:----------|:-----------:|
| `subm_flag = .false.` | Configuration of submerged inflow; Disable to activate surface inflow | *Disabled* by default|
| `flt_off_sw = .false.` | Configuration of floating offtake; Disable to allow users to specify the absolute elevation of groundwater withdrawal under `outl_elvs`  | *Disabled* by default |
| `seepage = .true.` | Configuration of seepage processing | *Disabled* by default|


#### Input Data
Model inputs include a time-series of hourly meteorological conditions, 3 daily inflows, and 7 daily outflows in the [specified format and column order](http://aed.see.uwa.edu.au/research/models/GLM/downloads/AED_GLM_v2_0b0_20141025.pdf) (Table 1).

**Table 1**  Input data for simulation of *Groose Dhuenn* Reservoir.

| Input Data     | Filename   | Code Section |
| ---------------- |:----------:|:-----------:|
| meteorological | *met_hourly_cc_kb_1996_2014.csv*| `&meteorology`|
| inflows | *inflow1_v2_1996-2015_moving_4d_meanKB.csv*, *inflow2_v2_1996-2015_moving_4d_meanKB.csv*, *inflow3_v2_1996-2015_moving_4d_meanKB.csv*,|`&inflow`|
| outflows | *outflowdep_1v2_1996-2014.csv*, *outflowdep_2v2_1996-2014.csv*, *outflowdep_3v2_1996-2014.csv*, *outflowdep_4v2_1996-2014.csv*, *outflowdep_5v2_1996-2014.csv*, *outflowdep_6v2_1996-2014.csv*, *outflow_loos_v2_1996-2015.csv*|`&outflow`|

### Example Output
Following each simulation, the model produces a .csv of time-series reservoir conditions detailing heat fluxes/water balance (*lake.csv*), an optional .csv of specific depth time-series water quality values (*WQ_2*), along with a NetCDF file detailing the values of each in-lake variables with time and depth (*output.nc*). For simple time-series plots, users can extract the values directly from the .csv file; For advanced plots, user may load the *output.nc* into MATLAB to produce the following outputs :

![alt text](https://eexlxq-dm2305.files.1drv.com/y4mr-ErYqZbV9smFf_eqiSCqLeJrHkUrLDnOrXtPHAfY8qGr_ZEZp_4Bymu0kJ6e6ltArqSuG3W25p5h-iVrutHi0VIuEkdUTMwEnihuyhRlY81aKdvEiEkMjGZZvAO3QVFftY8ZjX4Pjz-ddOlCVHmPzdQ2hS6CnrbgZuBtkIfkmOcSoJnotjg8TEbKl1nNmO6eWjBGoOTeIkA2pHEitw9vseqmJcyrAZx6DQkxDKA5Ug?width=1460&height=565&cropmode=none)
**Figure 1** Time-series Heat Fluxes of *Groose Dhuenn* Reservoir over 2 years

![alt text](https://fuxlxq.dm2301.livefilestore.com/y4m30FDfgm0jRh6l7GgpGCEAT-QL7KVtxPujU9pCWHzAbt8sjPVV7jgMYJ5wFruTB042qhWA9CzxCU1b8_F4LQKqoS2NYLgCyM-X2405cetrOAnsxcqvDGzaeT08SOcxIbDFTWGQXE0qS7yi3wRMCGgsbXnqk3ypDGtTdQuSF5x1q6xztkeZqigQ0UjH7j44PZp90sSPZwL_aKxTFOGLoYvFKcDA4GYvZQpHX3usM35xTI?width=1460&height=565&cropmode=none)
![alt text](https://foxlxq.dm1.livefilestore.com/y4mwXAiehDgxPb-WFmmwPlAMAudz_J4Nys2w4ur_7Yk8up1K2TEzpBloy67Wbgd0BvhUnhqI9Z5Qy6XfC5NSBUz0_A7m-RvyPe6IjrIBJr7rD0ndaLam-vd10aETzL7AMkH1LRlxTkG6FMr3myGrgbs2LhmtnNb-fMiT_mzoLl4iPuAxM_rRdTn9NEbfhJoh6PAgf5Bfk5U4iG1H_GU4oiS21iTtaS7dhaYDLfMLpm21tY?width=1460&height=565&cropmode=none)
**Figure 2** Time-series Water Balance of *Groose Dhuenn* Reservoir over 2 years


![alt text](https://epwaog-dm2305.files.1drv.com/y4mv9lo6ag_0XT6bfO2ttUR2xAIp_ODmwr-q4sOwLJNKHzpByaF1xjvNww06UpkGYSGA3VAqivj10bCHkSR-e5zze5YIVE-bCIi8-1xDFh7JLXikc3rFByLC25tze0sk71L2GMu8bZK7ysmlC8IaDjzUE7qM7CXdiKFK8SzzQ8dKqdKwOqgBIQEICQJxYBrqueAaXQM39EhdiUURSV_ijS7o4Dvj5DSV4MsNKZql_VvnAY?width=875&height=656&cropmode=none)
**Figure 3** Temperature-Depth contour plot of *Groose Dhuenn* Reservoir over 2 years


## Further Notes
- Users can switch off individual model components of interest by inserting `!` before the codes.
- Users can tweak the resolution of model simulation by adjusting `min_layer_thick` and `max_layer_thick`.
- Users can set output directory/filenames under `&output`, along with the option of exporting modelled data at different specified depths.
- Users can choose to incorporate AED2 water quality module by switching on `&wq_setup` component (i.e. remove `!` before the codes).
