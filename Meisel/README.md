# README
The Jupyter notebooks were made before the write up, so the most up to date information will be from the write up. Also I generally had to tweak my code when I transfered from a notebook to the main file in VS CODE. I was unable to finish the model.
## Here is a quick overview of what each file is
 - Charge_Transfer_init is a file that initializes all of my charge transfer values. All values were origionally initializes in the params document but it was getting too crowded so I moved the charge transfer values over to a new file
 - Electron_proton_oxygen_ion_conduction does pretty much what the title says it will do as it hashes out ionic and electrionic conduction
 - Same with Gas_diffusion_modeling
 - PCEC_Calculations contains miscilanius calculations to de-clutter the params file
 - PCEC_Charge_Transfer_model: Im pretty sure the good file is the one without any spaces (the second one).  This was my first notebook and hashes out a lot of the project.  The write up is a more consise version of the notes in this notebook
 - PCEC_projec: is the project. This contians the model function.
 - PCEC_params: the biggest file. This creates teh solution vector and the pars class. There are many parameters that are initialized here and there are also many calculations done here. This file could be organized better.
