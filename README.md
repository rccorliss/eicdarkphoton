# eicdarkphoton

This is the code to generate and interact with unsmeared and smeared
EIC data.


my current workflow for interacting with both the EIC smearing code
and my independent analysis involves three terminal windows:

#1) eic terminal: Runs the eic code and environment.
#Has to mount cvmfs if I have rebooted since last use
sudo mount -t cvmfs eic.opensciencegrid.org /cvmfs/eic.opensciencegrid.org ; sudo mkdir -p /cvmfs/eic.opensciencegrid.org
#then jump into vagrant and the singularity:
cd ~/singularity_vm # where the vagrant file lives.
vagrant up #starts this vm if it isn't up, returns a mild rebuke if it is.
vagrant ssh -- -X #(or just vagrant ssh.  The extra stuff is maybe necessary for xforwarding)
#now you are inside the vm.  
cd /osx_sphenix/Singularity #this is from inside the vm, and takes us to the mount where the .simg lives
singularity shell --writable --bind /osx_sphenix -B /cvmfs:/cvmfs /cvmfs/eic.opensciencegrid.org/singularity/rhic_sl7_ext  # start up the container, telling it to pass the /osx_sphenix and cvmfs mounts in from vagrant.
#now you are inside the singularity
source /cvmfs/eic.opensciencegrid.org/x8664_sl7/MCEG/releases/etc/eic_bash.sh # load EIC environmental variables
alias ll='ls -Al' #because I use this habitually. 
#now the stuff is set up the way I like it
cd eic_data_passthrough/
#and now we can run the code
root -l make_tree.C\(\"../sum100_eic20x250_ep_epee_m5GeV_th_1deglab.djangoh.txt\"\)
# or, rather, the shell script that parses things for me:

#2) passthrough:  handles the copying from the passthrough folder to
#the local folder, since there is an issue mounting the external folder.
#first jump into vagrant:
cd ~/singularity_vm # where the vagrant file lives.
vagrant up #starts this vm if it isn't up, returns a mild rebuke if it is.
vagrant ssh -- -X #(or just vagrant ssh.  The extra stuff is maybe necessary for xforwarding)
#now you are inside the vm and can do the copying back and forth, eg:
cp /osx_sphenix/eic_data_passthrough/*sh ./eic_data_passthrough/
cp /osx_sphenix/eic_data_passthrough/*.C ./eic_data_passthrough/
cp /osx_sphenix/eic_data_passthrough/*.djangoh* ./eic_data_passthrough/

3) local machine

ConvertAllMG.sh  [reldir]:  Convert all the madgraph .ttree.root files in the
given RELATIVE directory into .otree.root files and their associated
djangoh files.

Once they're made, they need to be moved into a
singularity by hand so that they can be parsed into smeared files and
brought back to my local machine.

ReadMGsimple.C:  Reads the .ttree.root files from madgraph.  This is
the code that ConvertAllMG.sh calls to convert files in the chosen directory.

GenerateReach.C:  Reads my (yeah, yeah, another format) 'otree' which
contains useful variables for A' reconstruction and generates various
plots, as well as reach calculations.

