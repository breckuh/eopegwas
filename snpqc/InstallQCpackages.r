###################################################################################################
#################### installing required packages to run snpQC ####################################
#################################################################################################

installIfNeeded = function (packages, installFn = install.packages) {
  newPackages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(newPackages)) installFn(newPackages)
}

installIfNeeded(c("gplots","gtools","lattice","RSQLite","snowfall","gdata"))



# snpQC uses a few packages - you can run the following line to install them:
install.packages(c("gplots","gtools","lattice","RSQLite","snowfall","gdata"),dependencies=T)
# and then choose a mirror to download from (a list should pop up)

# if you are running R without a graphical interface, use: 
#setRepositories(graphics=F, ind=1:2) # set repositories for CRAN
#local({r = getOption("repos"); r["CRAN"] = "http://cran.ms.unimelb.edu.au"; options(repos=r)}) # set to Melbourne mirror (or any mirror you prefer)
#install.packages(c("gplots","gtools","lattice","RSQLite","snowfall","gdata"),dependencies=T)

# if you are running on a server you might have problems with permissions to install packages - create a personal library, install the packages there and set the path to it 
# to install:
#myRlib="/home/user/Rlib" # for example
#install.packages(c("gplots","gtools","lattice","RSQLite","snowfall","gdata"),dependencies=T,lib="myRlib")

# and to use the packages, set the path with 
#.libPaths("/home/user/Rlib/") # before running snpQC
