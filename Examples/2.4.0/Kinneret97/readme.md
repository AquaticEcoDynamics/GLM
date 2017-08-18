# Lake Kinneret

## Study Site
Lake Kinneret is a [freshwater lake](https://www.researchgate.net/publication/289411199_Water_quality_modelling_of_Lake_Kinneret_Sea_of_Galilee) (32 ºN, 35 ºE) situated in northern Israel, and is the lowest freshwater lake in the world (209m below Mediterranean Sea Level).
The lake has a surface area of 173 km<sup>2</sup>, and is about 49 m deep.

## Model Setup
The model is set-up to simulate the hydrological domain of Lake Kinneret for 2 years, from **1997-01-01** to **1999-01-01** (`&time`), with water balance and heat fluxes [hypothetically calculated](http://aed.see.uwa.edu.au/research/models/GLM/downloads/AED_GLM_v2_0b0_20141025.pdf) based on the lake configuration and input data. Water quality functionality (AED2) is disabled for the simulation (`&wq_setup`).


#### Lake Configuration
The flexible model structure allow users to configure the switches of the individual model components, enabling the customization of the physical model, meteorological and flow conceptualisations.


**Physical Model Configuration**

| Configuration Flags | Description | Activation |
| ---------------- |:----------|:-----------:|
| `&morphometry` | Lake location (*latitude and longitude*) and bathymetry (*depth and area*) setup | - |
| `&wq_setup` | Water Quality setup | *Disabled* |
| `&init_profiles` | Initial lake conditions with depths (eg. *temperature and salinity* ) | - |
| `&output` | Specification of output file details | - |

**Meteorological Configuration**

| Configuration Flags | Description | Activation|
| ---------------- |:----------|:-----------:|
| `met_sw = .true.` | Surface meteorological forcing | *Enabled* |
| `subdaily = .true.` | Subdaily meteorological input setup switch (in *hourly timesteps*; `dt = 3600`) | *Enabled* |
| `&snowice` | Snow-ice lake dynamics setup | *Enabled* |
| `&bird_model` | Incorporation of *Bird Solar Radiation* model which theoretically computes solar radiation based on available data| *Enabled* |

**Flow Configuration**

| Configuration Flags | Description | Activation |
| ---------------- |:----------|:-----------:|
| `subm_flag = .false.` | Configuration of submerged inflow; Disable to activate surface inflow | *Disabled* by default|
| `flt_off_sw = .false.` | Configuration of floating offtake; Disable to allow users to specify the absolute elevation of groundwater withdrawal under `outl_elvs`  | *Disabled* by default |
| `seepage = .true.` | Configuration of seepage processing | *Disabled* by default|


#### Input Data
Model inputs include a time-series of hourly meteorological conditions, 6 daily inflows, and 1 daily outflow in the [specified format and column order](http://aed.see.uwa.edu.au/research/models/GLM/downloads/AED_GLM_v2_0b0_20141025.pdf) (Table 1).

**Table 1**  Input data for simulation of Lake Kinneret.

| Input Data     | Filename   | Code Section |
| ---------------- |:----------:|:-----------:|
| meteorological | *met_1997_2004.csv*| `&meteorology`|
| inflows | **Jordan** *inflow_1_v2_wq_1997_2004.csv*, **Neger** *inflow_2_v2_wq_1997_2004.csv*, **Yarmuch** *inflow_3_v2_wq_1997_2004.csv*, **Groundwater** *inflow_4_v2_wq_1997_2004.csv*, **Kinneret7** *inflow_5_v2_wq_1997_2004.csv*, **Salt Canal** *inflow_6_v2_wq_1997_2004.csv*, |`&inflow`|
| outflows | *outflow_v2_1997_2004.csv*|`&outflow`|

### Example Output
Following each simulation, the model produces a .csv of time-series reservoir conditions detailing heat fluxes/water balance (*lake.csv*), an optional .csv of specific depth time-series water quality values (*WQ_2*), along with a NetCDF file detailing the values of each in-lake variables with time and depth (*output.nc*). For simple time-series plots, users can extract the values directly from the .csv file; For advanced plots, user may load the *output.nc* into MATLAB to produce the following outputs :

![alt text](https://fuybrg.dm2302.livefilestore.com/y4mrPpzJaxiaCzvSwj5P68EVVIcf3yqTTVOaFH711W4Wl51nyKPBk-2-X7YuRBbsHNo1d_TsHSqQ8ym2vYhdAXK1rghFA6943zWjUstpYTg7Vx1BtBuwSMFm66jrpKuoxL_GL0AIWfJh9m2J14pn0-Ix_OU3nYsT3RdcCbIfHtdGn-S2scpoWZ5EvVs03O9VLpiV9Fii-Gc6YD7NnuKpBw6h__15on4SaM0J9bohuaouAQ?width=1460&height=565&cropmode=none)
**Figure 1** Time-series Heat Fluxes of Lake Kinneret over 2 years

![alt text](https://gpybrg-dm2305.files.1drv.com/y4mHMR-dSFqcCnpBsJeDaDouiA_GYuT3wjrYkcTNVfeeENwCWRISAf9WASzFkvSWhdfIjdfWFDbVSw4972OJwdPHsShCBPb_fmyuO9vR8h-9IPHWp1RxD2TXHAq3-Z8TkM4EMoXI-9l8e7nznN4CyWHIvwgduC2JZTvbXFK6-8SvZtI242yKzZXphE2IcUq__4_lyKG-jFeid4DrkY0OxjJBn0TGj4jp-AoiW1WhApG5UI?width=1460&height=565&cropmode=none)
![alt text](https://guybrg-dm2306.files.1drv.com/y4m4HbqTyp2Qiuc2isUrvewt6LOXUb_xVkdYZepbfvwXH0L1YAM5lk7YWDqylsXQO1PFZGL8lAFZZRdLgbsuVJWAevfqA9s3P9TCBWaFAKh5p5MkWOG6_A1qaNA4NTF_RuBZWMRi2_E6giTyKoEP1if0ERmPYP75zLG6kAbSmAEiHJxDUNSg6O-AOaBgm7Ptaqx8WjL8eDolua_aganAyLMggBaFJt68Xqgy-y0fjnPsGk?width=1460&height=565&cropmode=none)
**Figure 2** Time-series Water Balance of Lake Kinneret over 2 years


![alt text](https://eow2iw-dm2pap001.files.1drv.com/y4mkDHLVtOzvbrZLIEPrav7tRuRcv0KfiBBB56pTGw7dAlPi-xr8lS9QatH4eTfcDTYHbloqGI0wHw1ZC6006BzO4yBiFp-GdhgewFZR6nxbg2y4ydpJVIOy0KQvrGbdoCO6-X4zbQprkgXDxgIQbyjFQXtqDE70zMNvIb90aNXQ5wu-y8ktdQZWmRzUw8YWqJO5-24SPyE_CBskXlqiQcg0QDJnXPOW0HfVIM14kjB0fQ?width=875&height=656&cropmode=none)
**Figure 3** Temperature-Depth contour plot of Lake Kinneret over 2 years


## Further Notes
- Users can switch off individual model components of interest by inserting `!` before the codes.
- Users can tweak the resolution of model simulation by adjusting `min_layer_thick` and `max_layer_thick`.
- Users can set output directory/filenames under `&output`, along with the option of exporting modelled data at different specified depths.
- Users can choose to incorporate AED2 water quality module by switching on `&wq_setup` component (i.e. remove `!` before the codes).
