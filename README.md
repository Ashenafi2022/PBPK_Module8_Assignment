# PBPK_Module8_Assignment


---
# title: "Module8_Assignment"
# author: "Ashenafi Beyi"
# date: "3/20/2021"
---
  install.packages("deSolve")
  
  library(deSolve)
# Physiological parameters
# Blood flow rates
QCC = 12.9  # Cardiac output (L/h/kg) (Brown et al., 1997, Table 22 for mixed-breed (Mongrel) dogs (8.39 L/h/kg) and beagles (12.9 L/h/kg))
QLC = 0.297 # Fraction of blood flow to the liver (Brown et al., 1997, Table 26)
QKC = 0.173 # Fraction of blood flow to the kidneys (Brown et al., 1997, Table 26)
QFC = 0.097 # Fraction of blood flow to the fat (Vinegar, 2001, Table 2)
QMC = 0.217 # Fraction of blood flow to the muscle (Brown et al., 1997, Table 26)

# Tissue volumes
BW = 11.3      # Body weight (kg) (Baggotet al., 1977)
VLC = 0.0329   # Fractional liver tissue (Brown et al., 1997, Table 6)
VKC = 0.0055   # Fractional kidney tissue (Brown et al., 1997, Table 6)
VFC = 0.15     # Fractional fat tissue (Vinegar, 2001, Table 2)
VMC = 0.4565   # Fractional muscle tissue (Brown et al., 1997, Table 6)
VbloodC = 0.082 #   Blood volume, fraction of BW (Brown et al., 1997, Table 21)

# Fraction of tissue volumes that is blood ! (Brown et al., 1997; Table 30)
FVBF = 0.02 # Blood volume fraction of fat (%)
FVBM = 0.01 # Blood volume fraction of muscle (%)
FVBS = 0.01 # Blood volume fraction of slowly perfused tissues (%)

# Mass Transfer Parameters (Chemical-specific parameters)
# Partition coefficients (PC, tissue:plasma)
PL = 1.89  # Liver:plasmaPC (Craigmillet al., 2000, Table 3; Craigmill, 2003, Table 4, in sheep)
PK = 4.75  # Kidney:plasmaPC (Craigmillet al., 2000, Table 3; Craigmill, 2003, Table 4, in sheep)
PM = 0.85  # Muscle:plasmaPC (Craigmillet al., 2000, Table 3; Craigmill, 2003, Table 4, in sheep)
PF = 0.086 #  Fat:plasmaPC (Craigmillet al., 2000, Table 3; Craigmill, 2003, Table 4, in sheep)
PR = 4.75  # Richly perfused tissues:plasmaPC (Assumed the same as kidney:plasmaPC)
PS = 0.85  # Slowly perfused tissues:plasmaPC (Assumed the same as muscle:plasmaPC)

# Permeability cinstant (L/h/kg tissue) (permeation area cross products)
PAFC = 0.012 # Fat tissue permeability (0.12 in Leavens et al., 2012)
PAMC = 0.225 # Muscle tissue permeability (0.45 in Leavens et al., 2012)
PASC = 0.049 # Slowly perfused tissue permeability (0.49 in Leavens et al., 2012)

# Kinetic constants
# Oral absorption rate constants
Kst = 2    #  /h, gastric emptying rate constant
Ka = 0.012 #  0.05 ; /h, intestinal absorption rate constant, ka=0.05 for experimental solution; ka=0.012 for tablets or capsules
Kint = 0.2 # /h, intestinal transit rate constant

# IM absorption rate constants
Kim = 0.3   #  0.15 for conventional formulation; 0.3 for long-acting formulation; IM absorption rate constant (/h)
Frac= 0.5   # 0.95 for conventional formulation; 0.5 for long-acting formulation
Kdiss= 0.02 # /h

# IV infusion rate constants
Timeiv= 0.01   # IV infusion time (h), based on Leavens et al., 2012

# Urinary elimination rate constant
KurineC= 0.2     # L/h/kg

# Urinary elimination rate
Kurine = KurineC * BW # L/h

# Parameters for various exposure scenarios
PDOSEoral= 100  # (mg/kg)
PDOSEiv= 0      # (mg/kg)
PDOSEim= 0      # (mg/kg)

# Cardiac output and blood flows to tissues (L/h)
QC = QCC*BW    # Cardiac output
QL = QLC*QC    # Liver
QK = QKC*QC    # Kidney
QF = QFC*QC    # Fat
QM = QMC*QC    # Muscle
QR = 0.626*QC-QK-QL # Richly perfused tissues
QS = 0.374*QC-QF-QM # Slowly perfused tissues

# Tissue volumes (L)
VL = VLC*BW    # Liver
VK = VKC*BW    # Kidney
VF = VFC*BW    # Fat
VM = VMC*BW    # Muscle
Vblood= VbloodC*BW  # Blood
VR = 0.142*BW-VL-VK # Richly perfused tissues
VS = 0.776*BW-VF-VM # Slowly perfused tissues

# Permeability surface area coefficients
PAF = PAFC*VF   # Fat:bloodpermeability (L/h)
PAM = PAMC*VM   # Muscle:bloodpermeability (L/h)
PAS = PASC*VS   # Slowly perfused tissue:bloodpermeability (L/h)

# Volume of tissue vs blood
# Fat
VFB= FVBF*VF  # Fat compartment blood volume
VFT= VF-VFB   # Fat compartment tissue volume

# Muscle
VMB = FVBM*VM  # Muscle compartment blood volume
VMT = VM-VMB   # Muscle compartment tissue volume

# Slowly perfused tissue
VSB = FVBS*VS  # Slowly perfused compartment blood volume
VST = VS-VSB   # Slowly perfused compartment tissue volume

# Dosing
DOSEoral = PDOSEoral*BW  # (mg)
DOSEiv = PDOSEiv*BW      # (mg)
DOSEim = PDOSEim*BW      # (mg)

# Dosing, multiple oral gavage
DOSEimfast = DOSEim * Frac # Dose allocated to the fast absorption phase (mg)
DOSEimslow = DOSEim * (1 - Frac) # Dose allocated to the depot for slow absorption (mg)

# Dosing, multiple oral gavage
tlen= 0.1      # Length of oral gavage exposure (h/day)
tinterval= 6   # Varied dependent on the exposure paradigm
tdose = 1      # Number of oral gavage
OralR = DOSEoral/tlen

# OTC iv injection to the venous
IVR = DOSEiv/Timeiv

pbpkmodel <- function(Time, State, Parameters){
  with(as.list(c(State, Parameters)), {
    #########################################################################
    ## Concentration of the chemical in the vein of each compartment
    CVL = AL/(VL * PL) # Concentration in the liver venous blood
    CVK = AK/(VK * PK) # Concentration in the kidney venous blood
    CVR = AR/(VR * PR) # Concentration in the venous blood of richly perfused compartment
    CVF = AFB/VFB      # Concentration in the venous blood of fat
    CVS = ASB/VSB      # Concentration in the venous blood of slowly perfused compartment
    CVM = AMB/VMB      # Concentration in the venous blood of muscle
    
    ## Concentration of the chemical in the tissue sub-compartment of each member-limited compartment
    CMT = AMT/VMT  # Muscle
    CFT = AFT/VFT  # Fat
    CST = ASLT/VST # Slowly perfused tissues
    
    # OTC iv injection to the venous
    RIV = IVR * (Time < Timeiv)
    dAIV = RIV
    
    # OTC oral gavage to the stomach
    RDOSEoral <- OralR * (Time <= tdose * tinterval) * (Time %% tinterval < tlen) # Oral exposure rate (mg/h)
    dADOSEoral = RDOSEoral  # Amount that orally administered (mg)
    RAST = RDOSEoral - Kst * AST # Rate of change in the amiunt of chemical in the stomach (mg/h)
    dAST = RAST # Amount in the stomach (mg)
    
    ## OTC injection in intestine and colon
    RAI = Kst * AST - Kint * AI - Ka * AI # Rate of change in the amount of chemical in the intestine (mg/h)
    dAI = RAI # Amount in the intestine (mg)
    RColon = Kint * AI  # rate of change in the amount of chemical in the colon (mg/h)
    dAColon = RColon # Amount in the colon (mg)
    RAO = Ka * AI # Oral absorption rate (mg/h)
    dAAO = RAO # Amount absorbed orally (mg)
    
    # OTC in injection to the muscle
    RDOSEimremain = -Kdiss * DOSEimremain
    dDOSEimremain = RDOSEimremain
    Rim = Kim * Amtsite # im absorption rate (mg/h)
    dAbsorb = Rim # Amount absorbed after im injection (mg)
    Rsite = Rim + Kdiss * DOSEimremain # Rate of change in the amount of absorbable OTC in the injection site (mg/h)
    dAmtsite = Rsite # Amount of absorbable OTC that remains in the injection site (mg)
    
    #####################################################################################
    ## OTC in blood compartment
    CV = ((QL * CVL + QK * CVK + QF * CVF + QM * CVM + QR * CVR + QS * CVS + RIV + Rim)/QC) # Concentration in the vein
    CA = AA/Vblood # Concentration in the artery
    RA = QC * (CV - CA) # Rate of change in the amount of chemical in the blood compartment
    dAA = RA
    dAUCCV = CV
    
    ######################################################################################
    ## OTC in liver compartment
    RL = QL * (CA - CVL) + RAO # Rate of change in the amount of chemical in liver
    dAL = RL # Amount of chemical in liver
    CL = AL/VL
    dAUCCL = CL
    
    ################################################################################
    #OTC in kidney compartment
    # Urinary excretion of OTC
    Rurine = Kurine * CVK
    dAurine = Rurine
    
    # Kidney
    RK = QK * (CA - CVK) - Rurine # Rate of change in the amount of chemical in kidney
    dAK = RK # Amount of chemical in kidney
    CK = AK/VK # Concentration in kidney
    dAUCCK = CK
    
    #################################################################################
    ## OTC in muscle compartment, permeability-likited model
    RMB = QM * (CA - CVM) - PAM * CVM + (PAM * CMT)/PM # Rate of change in the amount of chemical in blood sub-compartment
    dAMB = RMB # amount of chemical in blood sub-compartment
    RMT = PAM * CVM - (PAM * CMT)/PM # Rate of change in the amount of chemical in tissue sub-compartment
    dAMT = RMT # amount of chemical in tissue sub-compartment
    AMtotal = AMT + AMB # total amount of chemical in muscle
    CM = AMtotal/VM # average/overall concentration of chemical in muscle
    dAUCCM = CM
    
    ##########################################################################
    ## OTC in fat compartment, permeability-limited model
    RFT = PAF * CVF - (PAF * CFT)/PF # Rate of change in the amount of chemical in tissue sub-compartment
    dAFT = RFT # Amount of chemical in tissue sub-compartment
    
    RFB = QF * (CA - CVF) - PAF * CVF + (PAF * CFT)/PF # Rate of change in the amount of chemical in blood sub-compartment
    dAFB = RFB # amount of chemical in blood sub-compartment
    AFtotal = AFT + AFB # total amount of chemical in fat
    CF = AFtotal/VF  # average/overall cncentration of chemical in fat
    
    ########################################################################
    ## OTC in richly perfused tissue compartment
    RR = QR * (CA - CVR) # Rate of change in the amount of chemical in richly perfused tissues
    dAR = RR  # Amount of chemical in richly perfused tissues
    CR = AR/VR # Concentration in richly perfused tissues
    
    ########################################################################
    ##  OTC in slowly perfused tissue compartment, permeability-limited nodel
    RSB = QS * (CA - CVS) - PAS * CVS + (PAS * CST)/PS # Rate of change in the amount of chemical in blood sub-compartment
    dASB = RSB # amount of chemical in blood sub-compartment
    RSLT = PAS * CVS - (PAS * CST)/PS # Rate of change in the amount of chemical in tissue sub-compartment
    dASLT = RSLT # amount of chemical in tissue sub-compartment
    AStotal = ASLT+ ASB # total amount of chemical in slowly perfused tissues
    CS = AStotal/VS # average/overall concentration of chemical in slowly perfused tissues
    
    #########################################################################
    ## Mass balance
    Qbal = QC - QL - QK - QM - QF - QR - QS
    Tmass = AA + AL + AK + Aurine + AMtotal + AFtotal + AR + AStotal
    list(c(dAbsorb, dAmtsite, dDOSEimremain, dAIV, dADOSEoral, dAST,
           dAI, dAColon, dAAO, dAA, dAL, dAK, dAurine, dAMB,
           dAMT, dAUCCV, dAUCCL, dAUCCK, dAUCCM, dAFB, dAFT,
           dAR, dASB, dASLT))
  })
}

State <- c(Absorb = 0, Amtsite = DOSEimfast, DOSEimremain = DOSEimslow, AIV = 0, ADOSEoral = 0, AST = 0,
           AI = 0, AColon = 0, AAO = 0, AA = 0, AL = 0, AK = 0, Aurine = 0, AMB = 0,
           AMT = 0, AUCCV = 0, AUCCL = 0, AUCCK = 0, AUCCM = 0, AFB = 0, AFT = 0,
           AR = 0, ASB = 0, ASLT = 0)
Parameters <- c(Kst, Ka, Kint, Kim, Kdiss, QC, QL, QK, QM,QF, QR, QS, PAF, PAM, PAS)

#Time parameters
starttime <- 0
stoptime <- 12
dtout <- 0.1
Times <- seq(starttime, stoptime, dtout)

# PBPK OUTPUT
OUT <- ode(y = State,
           times = Times,
           func = pbpkmodel,
           parms = Parameters,
           method = 'lsoda')
pbpkout <- as.data.frame(OUT)

# Check mass balance
Tmass = pbpkout$AA + pbpkout$Al + pbpkout$AK + pbpkout$Aurine + pbpkout$AMT + pbpkout$AMB +
  pbpkout$AFT + pbpkout$AFB + pbpkout$AR + pbpkout$ASLT + pbpkout$ASB 
Bal = pbpkout$AAO + pbpkout$AIV + pbpkout$Absorb - Tmass
plot(Times, Bal, col="red")



Additonal plot



Dstart= 0.0    # Initiation day of oral gavage (day)
Dstop= 0.2     # Termination day of oral gavage (day)
MAXT = 1.0     # maximum comm. interval
CINTC = 0.1    # Communication interval
CINT = CINTC   # Communication interval
Tsim= STOPTIME #  Tstopin hours
DS = Dstart*24 # Initiation time point of oral gavage (h)
Doff = (Dstop-Dstart)*24 # Oral gavage duration (h)
TimeOn= Dstart*24        # Initiation time point of oral gavage (h)
TimeOff= Dstop*24+tlen   # Termination time point of oral gavage (h)

# Exposure = PULSE(0,tinterval,tlen)*PULSE(DS, Tsim, Doff)
Exposure1 = SQUAREPULSE(0,tlen)
Exposure2 = SQUAREPULSE(0,tlen) + SQUAREPULSE(0+tinterval,tlen)
Exposure5 = SQUAREPULSE(0,tlen) + SQUAREPULSE(0+tinterval,tlen) + SQUAREPULSE(0+2*tinterval,tlen) + SQUAREPULSE(0+3*tinterval,tlen) + SQUAREPULSE(0+4*tinterval,tlen)
