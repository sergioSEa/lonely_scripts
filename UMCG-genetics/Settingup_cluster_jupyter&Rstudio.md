Some instructions to set up *jupyter*  
`conda activate "YOUR_ENV_NAME"`  
`conda install -c anaconda jupyter`  

`jupyter notebook --generate-config`  
`vim  ~/<your_username>/.jupyter/jupyter_notebook_config.py`  

Uncomment these lines from jupyter_notebook_config.py  and modify them (or write them if not present  
> c.NotebookApp.allow_origin = '*'  
> c.NotebookApp.ip = '0.0.0.0'  

Set up a password  
`jupyter notebook password`
Start jupyter
`jupyter notebook  --ip=0.0.0.0 --no-browser --port “HIGH_NUMER: e.g. 8080”`  

Some instructions to set up *RStudio*  
These instructions require Singularity.  

In the folder where you want to keep the rstudio singularity container  
`mkdir -p run var-lib-rstudio-server`

`printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf`  
2. Copy the image /shares/CIBIO-Storage/CM/scratch/users/davide.golzato/test_singularity/tidyverse_latest_update.sif  or one from here (https://rocker-project.org/images/) in the folder with the container. Mine has system dependencies already included that weren't present in the ones from the website and you should be able also to install packages from github.  


3. Build the container  
Note: To use singularity, it should be installed in your path. PATH with singularity: export PATH=/shares/CIBIO-Storage/CM/scratch/tools/singularity-ce-3.11.4/bin:$PATH  
`singularity build --sandbox rstudio_instance /shares/CIBIO-Storage/CM/scratch/users/davide.golzato/test_singularity/tidyverse_latest_update.sif`  

These sif files are the core of singularity. It contains the whole environment, and can be transfer from one computer to another. To generate one, you can initiate an instance, install things there and then reverse create the sif file.  


4. Start the container instance. You can bind whatever CM path you want in the container's mountpoint /mnt/. So far I've just tried with my user folder in scratch.  
For instance, I have mount two different paths in two folders in /mnt/ : Syntax: Location_real,:Location_intance  

`singularity instance start --writable --bind run:/run,var-lib-rstudio-server:/var/lib/rstudio-server,database.conf:/etc/rstudio/database.conf,/shares/CIBIO-Storage/CM/scratch/projects/sandreu_phylogenies/:/mnt/project,/shares/CIBIO-Storage/CM/scratch/users/sergio.andreusanchez/:/mnt/home rstudio_instance/ rstudio_instance`  


5.  Inside a screen, launch this command  
`singularity shell --writable rstudio_instance instance://rstudio_instance`  


6. Launch the actual rstudio server session  

`PASSWORD="pass123" /usr/lib/rstudio-server/bin/rserver --auth-none=0 --auth-pam-helper-path=pam-helper --server-user=$(whoami) --auth-timeout-minutes=0 --auth-stay-signed-in-days=30 --www-port 8786`   

The port you choose are the last numbers.  

7. Exit the screen and access from the browser at the address:  
`http://{cm3,cm4,cm5}.cibio.unitn.it:<open_port_you_chose>`

