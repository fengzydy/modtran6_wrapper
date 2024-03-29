Jing Feng (jing.feng3@mail.mcgill.ca, jing.feng@noaa.gov), Oct 10, 2022

# Step1: 

Download MODTRAN6 (http://modtran.spectral.com) and purchases license.

# Step2: 

For students and researchers at McGill AOS:

A float license server has been installed, there is no need to activate it in your own workstation.

Please contact Yi Huang (yi.huang@mcgill.ca) for permission to the license server.

Modify your shell environment .bashrc file for connection to float license server.

Restart your shell to activate the shell environment.
       
# Step3: 

Check 'run_modtran_example.m' for an example code to run modtran6.

IMPORTANT: Please read through MODTRAN documentation to make sure tape files are written properly.

Tips: Use tp5tojson to convert *.tp5 file to *.json file for checking input format and values being assigned to each card. 
      This *.json file can be used alternatively as an input file to MODTRAN6.
