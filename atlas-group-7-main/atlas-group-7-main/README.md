# Setting up the code 
## From gitlab repository:
This is the instruction to setup the git repository for the ATLAS-FP2 experiment.

The experiment has four parts:
1. event display -- proceed to directory 'Atlantis'
2. energy calibration -- go to directory 'ZeeFit'
3. W-mass measurement -- go to directory 'Wmass'
4. Higgs discovery -- go to directory 'HiggsSearch'

Each group must work in their own copies of the main code provided. The copies are provided to you by your tutor (please ask them for access rights!). 

### Clone the repository

Now you can make a clone of this forked project so that you can work locally in your machines. Click 'Clone or download' and then 'Use HTTPS'. Now copy the https link you see. 

Open a fresh terminal and create a directory for your experiment:

 ```bash
 mkdir ATLASExperiment
 cd ATLASExperiment
 ```

Now do:
 ```bash
 git clone <paste  the https link you copied here>
 ```
for example, git clone https://gitlab.com/fp2atlas2022/fp-2-atlas-group-0.git

Now you have a local version of the code in your respective machines. 

**You have successfully completed setting up the code in your local machines.**

# Accessing data
**Copy** the data from the data folder to your local repository in a folder named data (maybe you have to produce it with mkdir data).

# Running the code

The following instructions are for ZeeFit, Wmass and HiggsSearch. For event display, you can directly proceed to the  Atlantis directory.



## Docker

From the ATLASExperiment directory, do
 ```bash
sudo docker run -it --rm -p 8888:8888 -v $(pwd):/home/fp_student/work gitlab-registry.cern.ch/ssolomon/docker-fp2atlas/fp-image-test:latest
 ```
Make sure your ATLASExperiment directory contains the downloaded ntuples and you run docker in the same directory.Please start your docker from the ATLASExperiment folder. 

The docker image is made from https://gitlab.cern.ch/ssolomon/docker-fp2atlas.

Copy and paste to your browser the link http://127.0.0.1:8888/?token=.......

**You are ready to start the experiment now!**

Please look also in the README files in the different subfolders. 

**In Linux system,** if you see permission issues associated with connection to docker daemon socket at unix, try the above command with sudo as:
 ```bash
sudo docker run -it --rm -p 8888:8888 -v $(PWD):/home/fp_student/work gitlab-registry.cern.ch/ssolomon/docker-fp2atlas/fp-image-test:latest
 ```
If the url asks you for password or token, you can either choose a password of yours or copy and paste the token from the http://127.0.0.1:8888/.... link



