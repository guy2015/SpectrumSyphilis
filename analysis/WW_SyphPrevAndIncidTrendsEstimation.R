#devtools::install_github("guy2015/SoSSyphTrends", auth_token="<YOUR TOKEN>")

rm(list=ls()); kkk = gc(verbose=FALSE);

require(SoSSyphTrends) # Load the package

year_predict=1990:2021


namesCol = c("Country","ISO3","ISO3_letters","WHO_region","Data_type","Data_type_code","Sex","Year"
             ,"Diagnostic_test","DX_Code","N_positive","N_tested","Prevalence",
             "Weight_for_Spectrum_fitting", "WghtNSpectrum","Pop2012_15to49","PoptimesWeights",
             "In_Multi_country_analysis", "N_eligible_for_testing",
             "Age", "Denominator_survey_GARPR","Nbr_Pregnancies_incl_outside_ANC_Newman2013",
             "Location_Sites","NUMBER_of_sites","Regional","Diagnostic_Desciption","Source","Comment",
             "Row_num_IHME_GBD_db_Nov2016","Additional_row_numbers_IHME")


TFnamesCol = c("Country",	"ISO3",	"ISO3_letters",	"WHO region",	"Data type",	"Data type, code",
               "Sex",	"Year",	"Diagnostic test",	"DX_Code",	 "N positive", 	 "N tested", 	"Prevalence",
               "Weight for Spectrum fitting",	"WghtNSpectrum",	"Pop2012 15-49y M+F",
               "WghtSpect*Pop2012",	"In Multi-country analysis?", "N eligible for testing",
               "Age",	"Denominator survey (> GARPR)",
               "# Pregnancies incl. outside ANC (Newman-L et al. 2013)",	 "Location / Sites",
               "NUMBER of sites",	 "Regional", 	"Diagnostic Desciption",	"Source",
               "Comment",	"Row_num, IHME GBD db Nov.2016",	"Additional row numbers, IHME","DiagTestAdjusteFactor")

options(java.parameters = "-Xmx100048m")

require(xlsx)

# Data base file names
name.data.file_OldRes = "STIDatabase SPECTRUM_06Feb201820180404_zerAdj_0.01Combined_out.xlsx"
name.data.file = "STIDatabase SPECTRUM_2018April17.xlsx"
name.popu.file = "PopSize_15to49_9021.XLSX"


SyphData = read.xlsx(name.data.file,sheetName="SyphData",startRow=1,#endRow=1638,
                     colIndex=1:length(namesCol), check.names=FALSE)

# Prevalence adjustment
#SyphData$Prevalence = 1.1*SyphData$Prevalence
DiagnosticTest = read.xlsx(name.data.file,sheetName="DiagnosticTests",startRow=1,endRow=25, colIndex=1:14,header=TRUE)

colnames(DiagnosticTest)=c("STI","STI_code",	"Specimen",	"Specimen_code",	"Sex",	"Diagnostic test",	"DX_code",	"Sensitivity",
                           "Source_for_sensitivity",	"Specificity",	"Source_for_specificity",	"Adjustment_factor",
                           "Source_for_adjustment_factor")

DiagnosticTest = DiagnosticTest[16:21,]

high_risk_adj = 1.1

SyphData$Prevalence = sapply(1:length(SyphData$Prevalence),function(ii){
  res = NA;
  if(!is.na(SyphData$DX_Code[ii]))
  {
    adj = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==SyphData$DX_Code[ii])];
    res = SyphData$Prevalence[ii]*adj*high_risk_adj
  }
  res
} )

SyphData$DiagTestAdjusteFactor = sapply(1:length(SyphData$Prevalence),function(ii){
  res = NA;
  if(!is.na(SyphData$DX_Code[ii]))
  {
    res = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==SyphData$DX_Code[ii])];
  }
  res
} )

#End  Prevalence adjustment
ISO3 = read.xlsx(name.data.file,sheetName="ISO3",startRow=1, colIndex=1:3,header=FALSE)
colnames(ISO3)=c("Counrty","Numeric_ISO3","Alpha_ISO3")

#names(SyphData) = namesCol;
namesCol = c(namesCol,"DiagTestAdjusteFactor")
SyphData = SyphData[,is.element(names(SyphData),TFnamesCol)];
genv = environment();
ll = sapply(1:ncol(SyphData), function(jj) if(is.element(names(SyphData)[jj],TFnamesCol))
{
  iiind = which(TFnamesCol==names(SyphData)[jj])
  names(genv$SyphData)[jj]= namesCol[iiind]
  jj;
})

SyphData$WghtNSpectrum = SyphData$Weight_for_Spectrum_fitting

data.syphilis = SyphData;

PopSize = read.xlsx(name.popu.file,sheetName="Population aged 15-49",startRow=2, colIndex=c(3,8:39),header=FALSE)
colnames(PopSize) =  c("ISO3", paste("Year_", 1990+(0:(ncol(PopSize)-2)),sep=""))

#Syphilis
Sypdur_B = 2.42 # B
Sypdur_C = 4.13 # (C)

duration.syphilis=data.frame(t.zone=c("A","B","C"),duration=c(1.28,Sypdur_B,Sypdur_C))


# Loading functions
fCountryAnalysis_glob(Nboots=400, fname.data.file = name.data.file, Fmaxknots=2, FB_ProjMax=3, zerprev_adj=1/100) # Running the code
fCountryAnalysis_globN(Nboots=400, fname.data.file = name.data.file, Fmaxknots=2, FB_ProjMax=3, zerprev_adj=1/100,fname.data.file_OldRes=name.data.file_OldRes) # Running the code



