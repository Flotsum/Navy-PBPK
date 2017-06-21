#PROGRAM: PHYSIM -- Physiologically Based Pharmacokinetic Modeling Demo 

# Developed for ACSL/PC Level 10 -- July 15, 1993
# by Harvey Clewell (KS Crump Group, ICF Kaiser Int'l., Ruston, LA)
# and Mel Andersen (Health Effects Research Laboratory, USEPA, RTP, NC)
# Converted to R by Eric Hack -- March 2017

# Added DE for arterial blood

library(deSolve)


######################################################################
#### Define the model/system of ODEs
######################################################################

parms <- 1
dynamic <- function (times, state, parms) {
 with(as.list(c(state, parms)), {

#----DYNAMIC:  Things that change with time go in this function definition

#----Set up the dosing parameters   
#----IV = Intravenous infusion rate (mg/hr)
   IVZONE <- if (times > TINF) {0} else {1}     # 0 if times >= Length of IV infusion, 1 otherwise
   IV <- IVR*IVZONE                             # IVR= Intravenous infusion rate (mg/hr)
   
#----Inhalation
CIZONE <- as.double(CC | (times < TCHNG))       # Open or closed chamber (CC): 0 if CC=False and t>=TCHNG, 1 otherwise
   CI <- AI/VCH*CIZONE                          # Concentration in air (mg/L) = Amount in inhaled air/Net chamber volume (L)
   # Convert mass units to vomumetric gas units under std conditions. At 25 deg.C. the conversion fro ppm is 24450 mg/L
   # e.g https://ntp.niehs.nih.gov/iccvam/SuppDocs/FedDocs/OECD/OECD-GD39.pdf
   CP <- CI*24450./MW                           # Concentration (ppm)
   # Concentration (ppm)= Concentration in air (mg/L)*(conversion factor)/Molecular weight
   
#----Compute all concentrations at the top (need these for rate equations)
#----Slowyly perfused tissues
   CVS <- AS/(VS*PS)        # Concentration in venous blood leaving slowly-perfused tissue (mg/L)
   CS <- AS/VS              # Concentration in slowly-perfused tissues (mg/L)

#----Rapidly perfused tissues   
   CVR <- AR/(VR*PR)    # Concentration in venous blood leaving rapidly-perfused tissue (mg/L)
   CR <- AR/VR          # Concentration in rapidly-perfused tissue (mg/L)
   
#----Fat
   CVF <- AF/(VF*PF)    # Concentration in venous blood leaving fat tissue (mg/L)
   CF <- AF/VF          # Concentration in fat tissue (mg/L)

#----Liver
   CVL <- AL/(VL*PL)    # Concentration in venous blood leaving Liver (mg/L)
   CL <- AL/VL          # Concentration in Liver (mg/L)

#----CV = Mixed venous blood concentration (mg/L) 
   CV <- (QF*CVF + QL*CVL + QS*CVS + QR*CVR + IV)/QC 
   # Notice the IV dose entering the venous blood

#----CA = Concentration in arterial blood (mg/L)
   # CA <- (QC*CV+QP*CI)/(QC+(QP/PB))
   # Steady-state, rapid equilibration assumed (dAA/dt = 0)
   CA <- AA/VA                                  #Amount in Arterial blood/volume of arterial blood 
   
   CX <- CA/PB                                  # Concentration in exhaled air (mg/L)
   CXPPM <- (0.7*CX+0.3*CI)*24450./MW           # Concentration (ppm), 30% dead space assumed
   
#----AI = Amount in inhaled air (mg/L)
    dAI <- RATS*QP*(CA/PB-CI) - (KL*AI)         # Rate of change in air amount (mg/hr)
    dAInh <- QP*CI                              # Rate inhaled (mg/hr)
    
#----AX = Amount exhaled (mg) 
    dAX <- QP*CX                                # Rate equation (mg/hr)
    
#----AA = Amount in arterial blood (mg)
    dAA <- QC*(CV - CA) + QP*(CI - CA/PB)       # Rate equation (mg/hr)

#----AS = Amount in slowly perfused tissues (mg)
    dAS <- QS*(CA-CVS)                          # Rate equation (mg/hr)
    
#----AR = Amount in rapidly perfused tissues (mg)
    dAR <- QR*(CA-CVR)                          # Rate equation (mg/hr)
    
#----AF = Amount in fat tissue (mg)
    dAF <- QF*(CA-CVF)                          # Rate equation (mg/hr)
    
#----MR = Amount remaining in stomach (mg)
#    dMR <- -KA*MR 
     MR <- DOSE*exp(-KA*times)	                # analytical solution
    
#----AO = Total mass absorbed from stomach (mg)
    dAO <- KA*MR        # used as liver input, but no integrated to get AO
                        # rate of absorption * amount remaining in stomach
     AO <- DOSE-MR
    
#----AM = Amount metabolized (mg) 
    dAM <- (VMAX*CVL)/(KM+CVL) + KF*CVL*VL 
    
#----AL = Amount in liver tissue (mg)
    dAL <- QL*(CA-CVL) - dAM + dAO 


    #----Mass balance equations
    #----TMASS = total mass stored and remaining in stomach (mg) 
    TMASS <- AF+AL+AS+AR+MR
    
    #----OUT = amount removed by metabolism and exhalation (mg)
    OUT <- AM + AX
    
    #----DOSEX = Amount absorbed (mg)
    DOSEX <- AInh+AO+IVR*TINF
    
    #----MASSBAL = mass balance = in - out - stored
    MASSBAL <- DOSEX - OUT - TMASS
    #print(MASSBAL)
    
#----model function must return a list for the ode solver
     list(c(dAI, dAInh, dAX, dAS, dAR, dAF, dAM, dAL, dAA))
     
})
}

#END       ! End of dynamic section


######################################################################
#### Default constant parameter values
######################################################################

#-------Chemical parameters
     MW <- 104.0  # Molecular weight (g/mole)
     
#-------Physiological parameters

    QPC <- 14.    # Alveolar ventilation rate (L/hr)
    QCC <- 14.    # Cardiac output (L/hr)
    QLC <- 0.25   # Fractional blood flow to liver
    QFC <- 0.09   # Fractional blood flow to fat
     BW <- 0.22   # Body weight (kg)
    VLC <- 0.04   # Fraction liver tissue
    VFC <- 0.07   # Fraction fat tissue
    VAC <- 0.4	  # Fraction arterial blood

    # Partition coefficients
    PL <- 3.46			# Liver
    PF <- 86.5			# Fat
    PS <- 1.16			# Slowly perfused
    PR <- 3.46			# Rapidly perfused
    PB <- 40.2			# Blood:air
    
    # Metabolism parameters
    VMAXC <- 8.4
    KM <- 0.36			# Half-max concentration (mg/L); not scaled by BW
    KFC <- 0.
    
#-------Dosing parameters

  PDOSE <- 0.     # Oral dose (mg/kg)
 IVDOSE <- 0.     # IV dose (mg/kg)
   TINF <- .01    # Length of IV infusion (hrs)
   CONC <- 1000.  # Inhaled concentration (ppm)
     CC <- TRUE  # Default to open chamber; if true, Source animal parameters or defaults
  NRATS <- 3.     # Number of rats (for closed chamber)
    KLC <- 0.     # First order loss from closed chamber (/hr)        
   VCHC <- 9.1    # Volume of closed chamber (L)

#-------Timing commands: change these to increase or decrease exposure

  TSTOP <- 24.   # Length of experiment (hrs)
  TCHNG <- 24.    # Length of inhalation exposure (hrs)
 POINTS <- 25    # Number of points in plot
  CINT  <- TSTOP/POINTS  # Communication interval (for output, hrs)


######################################################################
# *** Modify constants ***
######################################################################
  
  # Source a file with modifications to the default constant parameter values.
  # This needs to be done before calculations based on the constants.
  # For example, load a rat or human physiological parameter file,
  #  or change the dosing parameters, or change the chemical.

  # Can have multiple R scripts called here. 
  
  #BUT BE CAREFUL OF THE ORDER!
  #(that is, when changing the same parameter in multiple files)


  # Set species-specific parameters (Remove # and modify file name accordingly)
        source("rat")          # if inhalation, Check CC = TRUE
        #source("human.R")       # if inhalation, Check CC = FALSE
  
  # Set chemical-specific parameters  (Remove # and modify file name accordingly)
        source("Styrene.R") 
        #source("MeBr.R")


#####################################################################
# Scaled or Adjusted or Calculated parameters
#####################################################################
 
# Closed chamber simulation accounts for losses in closed chamber experiments, which tended to be
# more substantial in older studies.
  
if (CC) { # Closed chamber simulation 
	RATS <- NRATS    
	KL <- KLC
} else {
	RATS <- 0.  
	KL <- 0. # (Turn off chamber losses so concentration remains constant)
}
      
if (PDOSE == 0.0) KA <- 0.  # Parenteral dosing
   
VCH <- VCHC-RATS*BW        # Net chamber volume (L)
AI0 <- CONC*VCH*MW/24450.  # Initial amount in chamber (mg)


#------Allometric scaling
# Flows (blood and respiration)
QC <- QCC*BW**0.74		# Cardiac output (L/hr)
QP <- QPC*BW**0.74		# Alveolar ventilation (L/hr)
QL <- QLC*QC		        # Liver blood flow (L/hr)
QF <- QFC*QC		        # Fat blood flow (L/hr)
QS <- 0.24*QC-QF		# Slowly perfused blood flow (L/hr)
QR <- 0.76*QC-QL		# Rapidly perfused blood flow (L/hr)
## What are 0.24 and 0.76?

# Tissue Volumes
VL <- VLC*BW		# Liver (L)
VF <- VFC*BW		# Fat (L)
VA <- VAC*BW		# Arterial blood (L)
VS <- 0.82*BW-VF	# Slowly perfused (L)
VR <- 0.09*BW-VL	# Richly perfused (L)
## Why does this not add up to 1?

# Doses
DOSE <- PDOSE*BW		# Oral dose (mg)
IVR <- IVDOSE*BW/TINF	        # IV dose RATE (mg/hr)
KA <- 0.			# Oral absorption rate constant (/hr)

# Metabolism
VMAX <- VMAXC*BW**0.7	        # Maximal rate (mg/hr)
KF <- KFC/BW**0.3		# First order elimination (/hr)


#####################################################################
# Integration

# States (The initial (state) values for the ODE system)
state <- c(AI = AI0,
           AInh = 0.,
           AX = 0., 
           AS = 0., 
           AR = 0., 
           AF = 0., 
           AM = 0., 
           AL = 0.,
           AA = 0.)
#,
#           MASSBAL = 0)

# Integration times
times <- seq(0, TSTOP, by = CINT)

# Call the integrator
out <- ode(y = state, times = times, func = dynamic, parms = parms)
y <- as.data.frame(out)


#####################################################################
# Plotting routines

   source("load.Sty.data.R")
   source("plot.Sty.R")

plot(out)

#END       ! End of program 
