# ocean-co2-absorption

## Intro
The ocean has absorbed 41% of all anthropogenic carbon to date, and this sink is critical for climate change going forward. Monitoring the carbon sink is critical for supporting global-scale management of the climate via managing the atmospheric CO2 concentration, i.e. the UNFCCC Paris Agreement. Ocean data are quite sparse and CO2 in water cannot be directly measured from space, thus intend to use machine learning approaches informed by oceanographic understanding to extrapolate from sparse data to global spatially resolved fields at at least monthly temporal resolution.

## Data Cleaning and Preprocessing
* preprocess/

## Data Analysis and Visualization
* notebooks/reading_SOCAT_data_example.ipynb: Example of reading SOCAT data. This notebook downloads data from figshare, extracts the data with xarray, and plots it with cartopy.
* notebooks/pCO2_testbed_member_001_data_visualization.ipynb: Data visualization on pCO2 testbed data.

## Models
* model/landschutzer_SOM-FFN/: Landschutzer SOM-FFN neural network model in Matlab.