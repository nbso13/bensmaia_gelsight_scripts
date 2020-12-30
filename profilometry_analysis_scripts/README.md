# bensmaia_gelsight_scripts

Library:

Profilometry struct analysis functions:
- bruteCropFit gets the min length in each direction and crops input structs to that minimum length.
- characterizeFilterFull visualizes data filtered and unfiltered, conducts fourier analysis and produces frequency response. 
- checkSizeMatch takes in two gel structs and checks to make sure size matches, so they can go through charactFilter
- compare_lines plots an overlapping plot of the two cross sections (used to compare gain values
- cropProfile crops the profilometry profile and modifies struct.
- gainAnalysis takes in a gel and no gel of the same image and compares the ratio of gel to no gel in terms of height at a given resolution
- profilometry2shape transforms profilometry into a touchsim shape and pin offset
- resampleToMin resamples each profile down to the min resolution at each dimension
- sectionbyMm takes in two x,y coordinate pairs, one lower left hand one upper right hand, in mm, and a gel or no_gel profilometry struct. Returns the profile between them.
- shape2profilometry transforms touchsim shape and offset into profilometry struct
- skinModel takes in a touchsim stimulus and applies the skin mechanics model, giving back the skin profile for that stimulus.
- visualizeProfile displays given profilometry struct
- freqAnalysisFull takes the full fourier transform

Touchsim Functions:
- surfTouchSim takes in a shape and pin offset pair and plots a surface
- shape2profilometry transforms touchsim shape and offset into profilometry struct
- skinModel takes in a touchsim stimulus and applies the skin mechanics model, giving back the skin profile for that stimulus.
- profilometry2shape transforms profilometry into a touchsim shape and pin offset
- TouchSimSkin takes in two profilometry structs and returns the gel but downsampled as the rest of them are, the "new gel" or the touchsim version of the skin and the new_no_gel, downsampled to touchsim level.

Scripts:
- main_touchsim_pipeline : main script for touchsim analysis
- measure_gain_charles: measuring profilometer gain
- measuring_conformance: measure conformance using craig stimuli
- process_profilometry_data: process profilometry data from csv to matlab struct
- translating_from_touchsim_to_profilometry_and_back: converting between shape and profilometry