Code for the Master Research Project of Yannick Badoux for the Astronmy & Data Science Master at Leiden University.

Project working title:
Does a planet keep its moons when it is ejected from the system during an encounter with another star?

File descriptions:
- initial_conditions.py: file to generate initial conditions
- run_simulation.py: file used to evolve the system
- analyse_result.py: functions to analyse the evolved system (e.g. check which particles are bound to each other)
- HB83.py: script that reproduces the cross section calculations from Hut & Bahcall 1983
- acs_binary.py: improved cross section determination algorithm. Comparison to the results from HB83.
- moviemaker.py: can be used to create a video of the simulation if you set plot to True in run_simulation
- cross_section_analysis.ipynb: notebook to analyze and plot the output from HB83.py and acs_binary.py

Example shell scripts showing how to run a cross section study as a function of v_inf can be found in the shell_scripts folder.