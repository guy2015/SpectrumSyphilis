#Extra plot

#' Generate plots for results obtained from the function RunFitSyphilis
#'
#' @param syphfits list of class "Syph-fit" (containing a file name and the data used to run the fits).
#' @param countries vector containing the names of countries for which charts displaying estimates are to be drawn.
#' @return Creates and saves two pdf files in the working directory. The created pdf files have names with the prefix "R_Proj_Prev_", and display Syphilis prevalence data and estimated trends by sex and risk.
#' @examples Not available
SyphPlots <- function(syphfits,countries=NULL, yplot=2000:2023)
{
  library(tidyverse)
  library(scales)
  library(gridExtra)

  mout_file <- syphfits$filename;
  SyphData <- syphfits$SyphData

  all_res <- openxlsx::read.xlsx(mout_file,sheet="SYPH_RBootstrap_All")
  all_res$datatype <- "Model"
  all_res$weight <- 1

  #Prevalence rate
  names(all_res)[names(all_res)=="MedianPrevF"] = "PrevMed_Females"
  names(all_res)[names(all_res)=="EstimatePrevF"] = "PrevEst_Females"
  names(all_res)[names(all_res)=="PrevLB_2.5%F"] = "PrevLB_Females"
  names(all_res)[names(all_res)=="PrevUB_97.5%F"] = "PrevUB_Females"

  names(all_res)[names(all_res)=="MedianPrevM"] = "PrevMed_Males"
  names(all_res)[names(all_res)=="EstimatePrevM"] = "PrevEst_Males"
  names(all_res)[names(all_res)=="PrevLB_2.5%M"] = "PrevLB_Males"
  names(all_res)[names(all_res)=="PrevUB_97.5%M"] = "PrevUB_Males"

  names(all_res)[names(all_res)=="MedianPrevM+F"] = "PrevMed_BothSexes"
  names(all_res)[names(all_res)=="EstimatePrevM+F"] = "PrevEst_BothSexes"
  names(all_res)[names(all_res)=="PrevLB_2.5%M+F"] = "PrevLB_BothSexes"
  names(all_res)[names(all_res)=="PrevUB_97.5%M+F"] = "PrevUB_BothSexes"

  #Case prevalence
  names(all_res)[names(all_res)=="CasePrev_EstF"] = "CasePrevEst_Females"
  names(all_res)[names(all_res)=="CasePrev_MedF"] = "CasePrevMed_Females"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%F"] = "CasePrevLB_Females"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%F"] = "CasePrevUB_Females"

  names(all_res)[names(all_res)=="CasePrev_EstM"] = "CasePrevEst_Males"
  names(all_res)[names(all_res)=="CasePrev_MedM"] = "CasePrevMed_Males"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%M"] = "CasePrevLB_Males"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%M"] = "CasePrevUB_Males"

  names(all_res)[names(all_res)=="CasePrev_EstM+F"] = "CasePrevEst_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_MedM+F"] = "CasePrevMed_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%M+F"] = "CasePrevLB_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%M+F"] = "CasePrevUB_BothSexes"

  #Incidence rates
  names(all_res)[names(all_res)=="EstimateIncF"] = "InciMed_Females"
  names(all_res)[names(all_res)=="MedianIncF"] = "InciEst_Females"
  names(all_res)[names(all_res)=="IncLB_2.5%F"] = "InciLB_Females"
  names(all_res)[names(all_res)=="IncUB_97.5%F"] = "InciUB_Females"

  names(all_res)[names(all_res)=="EstimateIncM"] = "InciMed_Males"
  names(all_res)[names(all_res)=="MedianIncM"] = "InciEst_Males"
  names(all_res)[names(all_res)=="IncLB_2.5%M"] = "InciLB_Males"
  names(all_res)[names(all_res)=="IncUB_97.5%M"] = "InciUB_Males"

  names(all_res)[names(all_res)=="EstimateIncM+F"] = "InciMed_BothSexes"
  names(all_res)[names(all_res)=="MedianIncM+F"] = "InciEst_BothSexes"
  names(all_res)[names(all_res)=="IncLB_2.5%M+F"] = "InciLB_BothSexes"
  names(all_res)[names(all_res)=="IncUB_97.5%M+F"] = "InciUB_BothSexes"

  #Incidence cases
  names(all_res)[names(all_res)=="CaseInc_EstF"] = "CaseInciMed_Females"
  names(all_res)[names(all_res)=="CaseInc_MedF"] = "CaseInciEst_Females"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%F"] = "CaseInciLB_Females"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%F"] = "CaseInciUB_Females"

  names(all_res)[names(all_res)=="CaseInc_EstM"] = "CaseInciMed_Males"
  names(all_res)[names(all_res)=="CaseInc_MedM"] = "CaseInciEst_Males"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%M"] = "CaseInciLB_Males"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%M"] = "CaseInciUB_Males"

  names(all_res)[names(all_res)=="CaseInc_EstM+F"] = "CaseInciMed_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_MedM+F"] = "CaseInciEst_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%M+F"] = "CaseInciLB_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%M+F"] = "CaseInciUB_BothSexes"

  long_all_res <- data.frame();

  for(ctr in unique(all_res$Country))
  {
    if(!is.element(ctr,all_res$Country)) next
    temp_long_ctr <- subset(all_res,Country==ctr)

    long_ctr_men <- data.frame(Country = ctr,
                               sex = "Males",
                               datatype = "Model",
                               weight = 1,
                               Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                               indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                               Median = c(temp_long_ctr$CaseInciMed_Males,temp_long_ctr$InciMed_Males,temp_long_ctr$CasePrevMed_Males,temp_long_ctr$PrevMed_Males),
                               BestFit = c(temp_long_ctr$CaseInciEst_Males,temp_long_ctr$InciEst_Males,temp_long_ctr$CasePrevEst_Males,temp_long_ctr$PrevEst_Males),
                               Lower = c(temp_long_ctr$CaseInciLB_Males,temp_long_ctr$InciLB_Males,temp_long_ctr$CasePrevLB_Males,temp_long_ctr$PrevLB_Males),
                               Upper = c(temp_long_ctr$CaseInciUB_Males,temp_long_ctr$InciUB_Males,temp_long_ctr$CasePrevUB_Males,temp_long_ctr$PrevUB_Males)
    )


    long_ctr_women <- data.frame(Country = ctr,
                                 sex = "Females",
                                 datatype = "Model",
                                 weight = 1,
                                 Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                                 indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                                 Median = c(temp_long_ctr$CaseInciMed_Females,temp_long_ctr$InciMed_Females,temp_long_ctr$CasePrevMed_Females,temp_long_ctr$PrevMed_Females),
                                 BestFit = c(temp_long_ctr$CaseInciEst_Females,temp_long_ctr$InciEst_Females,temp_long_ctr$CasePrevEst_Females,temp_long_ctr$PrevEst_Females),
                                 Lower = c(temp_long_ctr$CaseInciLB_Females,temp_long_ctr$InciLB_Females,temp_long_ctr$CasePrevLB_Females,temp_long_ctr$PrevLB_Females),
                                 Upper = c(temp_long_ctr$CaseInciUB_Females,temp_long_ctr$InciUB_Females,temp_long_ctr$CasePrevUB_Females,temp_long_ctr$PrevUB_Females)
    )

    long_ctr_both <- data.frame(Country = ctr,
                                sex = "BothSexes",
                                datatype = "Model",
                                weight = 1,
                                Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                                indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                                Median = c(temp_long_ctr$CaseInciMed_BothSexes,temp_long_ctr$InciMed_BothSexes,temp_long_ctr$CasePrevMed_BothSexes,temp_long_ctr$PrevMed_BothSexes),
                                BestFit = c(temp_long_ctr$CaseInciEst_BothSexes,temp_long_ctr$InciEst_BothSexes,temp_long_ctr$CasePrevEst_BothSexes,temp_long_ctr$PrevEst_BothSexes),
                                Lower = c(temp_long_ctr$CaseInciLB_BothSexes,temp_long_ctr$InciLB_BothSexes,temp_long_ctr$CasePrevLB_BothSexes,temp_long_ctr$PrevLB_BothSexes),
                                Upper = c(temp_long_ctr$CaseInciUB_BothSexes,temp_long_ctr$InciUB_BothSexes,temp_long_ctr$CasePrevUB_BothSexes,temp_long_ctr$PrevUB_BothSexes)
    )

    long_ctr <- rbind(long_ctr_men,long_ctr_women,long_ctr_both)

    min_year <- min(long_ctr$Year)
    #Adding data points
    temp_ctr <- subset(SyphData,Country==ctr & Weight_for_Spectrum_fitting>0)
    ppui <- temp_ctr$Prevalence
    ppui[ppui<=0] = 1/100;

    temp_all_res <- data.frame(Country = ctr,
                               sex = NA,
                               datatype = temp_ctr$Data_type,
                               weight = temp_ctr$Weight_for_Spectrum_fitting,
                               Year=temp_ctr$Year,
                               indicator="PrevalenceRate",
                               Median = temp_ctr$Prevalence/100,
                               BestFit = temp_ctr$Prevalence/100,
                               Lower = NA,
                               Upper = NA
    )

    sdall <- sqrt(ppui/100*(1-ppui/100)/temp_ctr$N_tested)
    sdall[is.na(sdall)] <- 0
    temp_all_res$Lower <- pmax(temp_all_res$Median-1.96*sdall,0)
    temp_all_res$Upper <- pmin(temp_all_res$Median+1.96*sdall,1)

    temp_all_res$sex[temp_all_res$datatype=="ANC Routine screening"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="ANC Survey"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="BloodDonor Screening Men"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="FSW"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="MSM"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="Survey LowRisk Men"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="BloodDonor Screening Men + Women"] <- "BothSexes"
    temp_all_res$sex[temp_all_res$datatype=="Male Sex Workers"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="MSM + MSW combined"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="PWID-Female"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="PWID-Male"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="Survey LowRisk Women"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="Trans-Genders"] <- "BothSexes"
    temp_all_res$sex[temp_all_res$datatype=="BloodDonor Screening Women"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="Prisoners, Men"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="Prisoners, Women"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="Wives of PWID"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="Survey LowRisk Men+Women"] <- "BothSexes"

    ln_na <- unique(c(which(is.na(temp_all_res$BestFit)),which(is.na(temp_all_res$Median))))
    if(length(ln_na)>0)temp_all_res <- temp_all_res[-ln_na,]
    long_ctr <- rbind(long_ctr,subset(temp_all_res,Year>=min_year))

    long_ctr$datatype <- factor(long_ctr$datatype, levels=c("Model","ANC Routine screening","ANC Survey","BloodDonor Screening Men",
                                                            "FSW","MSM","Survey LowRisk Men","BloodDonor Screening Men + Women",
                                                            "Male Sex Workers","MSM + MSW combined","PWID-Female","PWID-Male",
                                                            "Survey LowRisk Women","Trans-Genders","BloodDonor Screening Women",
                                                            "Prisoners, Men","Prisoners, Women","Wives of PWID","Survey LowRisk Men+Women"))

    #Final merge:
    long_all_res <- rbind(long_all_res,long_ctr)
  }

  m_data <- long_all_res

  old_plot_aggregate <- function(df, mtitle)
  {
    df$Estimate <- df$BestFit
    df <- subset(df,Year%in%yplot)
    if(nrow(df)==0) return(NULL)

    df$datatype <- factor(df$datatype, levels=c("Model","ANC Routine screening","ANC Survey","BloodDonor Screening Men",
                                                "FSW","MSM","Survey LowRisk Men","BloodDonor Screening Men + Women",
                                                "Male Sex Workers","MSM + MSW combined","PWID-Female","PWID-Male",
                                                "Survey LowRisk Women","Trans-Genders","BloodDonor Screening Women",
                                                "Prisoners, Men","Prisoners, Women","Wives of PWID","Survey LowRisk Men+Women"))




    ii_x <- which(df$datatype=="Model" & df$sex!="BothSexes")
    df <- df[-ii_x,]
    print(
      df %>% filter(indicator=="PrevalenceRate")%>% ggplot(subset=.(Country==mtitle & datatype=="Model"), aes(x=Year, y=100*Estimate, ymin = 100*Lower, ymax= 100*Upper,color=datatype,fill=datatype))+
        geom_line(data =. %>% filter(Country==mtitle& datatype=="Model"), aes(Year, 100*Estimate)) +
        geom_ribbon(data =. %>% filter(Country==mtitle& datatype=="Model"),alpha=0.2,linetype=0)+
        geom_pointrange(data =. %>% filter(Country==mtitle& datatype!="Model"), position=position_jitter(w = 0.05, h = 0),#position=position_jitter(width=0.5),
                        linetype='solid', aes(size=weight),size=0.75)+
        expand_limits(y = 0) +
        labs(title = paste(mtitle," Syphilis prevalence trend among adults (15-49 y)", sep=","), x = "Year", y="Test-adjusted prevalence, active syphilis (%)") +
        theme_minimal() +
        theme(legend.position = "bottom")
    )
  }

  plot_aggregate <- function(df, mtitle)
  {
    df$Estimate <- df$BestFit
    df <- subset(df,Year%in%yplot)
    if(nrow(df)==0) return(NULL)

    df$datatype <- factor(df$datatype, levels=c("Model","ANC Routine screening","ANC Survey","BloodDonor Screening Men",
                                                "FSW","MSM","Survey LowRisk Men","BloodDonor Screening Men + Women",
                                                "Male Sex Workers","MSM + MSW combined","PWID-Female","PWID-Male",
                                                "Survey LowRisk Women","Trans-Genders","BloodDonor Screening Women",
                                                "Prisoners, Men","Prisoners, Women","Wives of PWID","Survey LowRisk Men+Women"))





    print(
      df %>% filter(indicator=="PrevalenceRate")%>% ggplot(subset=.(Country==mtitle & datatype=="Model"), aes(x=Year, y=100*Estimate, ymin = 100*Lower, ymax= 100*Upper,color=datatype,fill=datatype))+
        geom_line(data =. %>% filter(Country==mtitle& datatype=="Model"), aes(Year, 100*Estimate)) +
        geom_ribbon(data =. %>% filter(Country==mtitle& datatype=="Model"),alpha=0.2,linetype=0)+
        geom_pointrange(data =. %>% filter(Country==mtitle& datatype!="Model"), position=position_jitter(w = 0.05, h = 0),#position=position_jitter(width=0.5),
                        linetype='solid', aes(size=weight),size=0.75)+
        facet_wrap(~sex)+
        expand_limits(y = 0) +
        labs(title = paste(mtitle," Syphilis prevalence trend among adults (15-49 y)", sep=","), x = "Year", y="Test-adjusted prevalence, active syphilis (%)") +
        theme_minimal() +
        theme(legend.position = "bottom")
    )
  }

  by_country <- m_data %>%
    count(Country, BestFit, Median,Lower, Year,Upper, sex, datatype, indicator, wt = BestFit, name = "Estimate") %>%
    split(.$Country)

  out_file_new <- paste(paste("R_Proj_Prev_Both+M+F_",as.character(Sys.Date()), sep=""), ".pdf", sep="")
  suppressWarnings({
    pdf(out_file_new, h = 10, w = 18)
    Map(plot_aggregate, by_country, names(by_country))
    dev.off()
  }
  )

  out_file_old <- paste(paste("R_Proj_Prev_BothSexes_",as.character(Sys.Date()), sep=""), ".pdf", sep="")
  suppressWarnings({
    pdf(out_file_old, h = 10, w = 18)
    Map(old_plot_aggregate, by_country, names(by_country))
    dev.off()
  }
  )
}


#' Generate plots for results obtained from the function RunFitSyphilis
#'
#' @param name.data.file name of the file containing.
#' @param countries vector containing the names of countries for which charts displaying estimates are to be drawn.
#' @return Creates and saves two pdf files in the working directory. The created pdf files have names with the prefix "R_Proj_Prev_", and display Syphilis prevalence data and estimated trends by sex and risk.
#' @examples Not available
#' @keywords internal
#' @export
SyphPlotsFF <- function(name.data.file=fin_name.data.file,f_high_risk_adj = 1.1,in_year_predict=1990:2023)
{
  namesCol = c("Country","ISO3","ISO3_letters","WHO_region","Data_type","Data_type_code","Sex","Year"
               ,"Diagnostic_test","DX_Code","N_positive","N_tested","Prevalence",
               "Weight_for_Spectrum_fitting", "WghtNSpectrum")

  TFnamesCol = c("Country",	"ISO3",	"ISO3_letters",	"WHO region",	"Data type",	"Data type, code",
                 "Sex",	"Year",	"Diagnostic test",	"DX_Code",	 "N positive", 	 "N tested", 	"Prevalence",
                 "Weight_for_Spectrum_fitting",	"WghtNSpectrum")

  options(java.parameters = "-Xmx100048m")

  high_risk_adj = f_high_risk_adj;
  require(xlsx)

  # Data base file names
  #name.data.file_OldRes = "Database SPECTRUM_2018April17_forReview20180604_zerAdj_0.01Combined_out.xlsx" #If the analysis was already run, the data base edited and the result file is available, please modify this accordingly.
  #name.data.file = "STIDatabase SPECTRUM_2021Jan27.xlsx"
  #name.popu.file = "PopSize_15to49_9021.XLSX"

  SyphData = read.xlsx(name.data.file,sheetName="SyphData",startRow=1,endRow=3011,
                       colIndex=1:(1+length(namesCol)), check.names=FALSE)[,-c(7,9)]

  names(SyphData)[names(SyphData)=="Weight for Spectrum fitting"] <- "Weight_for_Spectrum_fitting"
  names(SyphData)[names(SyphData)=="Midpoint study year"] <- "Year"
  names(SyphData)[names(SyphData)=="Population type"] <- "Data type"
  names(SyphData)[names(SyphData)=="Population code"] <- "Data type, code"
  #Prevalence adjustment
  #SyphData$Prevalence = 1.1*SyphData$Prevalence
  DiagnosticTest = read.xlsx(name.data.file,sheetName="DiagnosticTests",startRow=1,endRow=25, colIndex=1:14,header=TRUE)

  colnames(DiagnosticTest)=c("STI","STI_code",	"Specimen",	"Specimen_code",	"Sex",	"Diagnostic test",	"DX_code",	"Sensitivity",
                             "Source_for_sensitivity",	"Specificity",	"Source_for_specificity",	"Adjustment_factor",
                             "Source_for_adjustment_factor")

  DiagnosticTest = DiagnosticTest[16:22,]

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

  #Syphilis
  Sypdur_B = 2.42 # B
  Sypdur_C = 4.13 # (C)

  temp_mout_file <- substr(mout_file,4,nchar(mout_file)-5)
  mout_file <- dir()[grep(temp_mout_file,dir())]
  if(length(mout_file)>=1)
  {
    mout_file <- mout_file[nchar(mout_file)==max(nchar(mout_file))]
    mout_file <- mout_file[length(mout_file)]
  }

  if(length(mout_file)==0)
  {
    mout_file <- dir()[grep("out",dir())]
    if(length(mout_file)>1) mout_file <- mout_file[length(mout_file)]
  }

  if(length(mout_file)==0) stop("No file available. Please check.")

  m_syphfits <- list(filename=mout_file,SyphData=SyphData)
  SyphPlots(m_syphfits)
}



#' Generate Country specific plots for results obtained from the function RunFitSyphilis
#'
#' @param syphfits list of class "Syph-fit" (containing a file name and the data used to run the fits).
#' @param countries vector containing the names of countries for which charts displaying estimates are to be drawn.
#' @return Creates and saves two pdf files in the working directory. The created pdf files have names with the prefix "R_Proj_Prev_", and display Syphilis prevalence data and estimated trends by sex and risk.
#' @examples Not available
np_SyphPlots <- function(syphfits,countries_iso=NULL, yplot=2000:2023)
{
  library(tidyverse)
  library(scales)
  library(gridExtra)

  #mout_file <- syphfits$filename;
  SyphData <- syphfits$SyphData

  all_res <- openxlsx::read.xlsx(syphfits$wb,sheet="SYPH_RBootstrap_All")
  all_res$datatype <- "Model"
  all_res$weight <- 1

  #Prevalence rate
  names(all_res)[names(all_res)=="MedianPrevF"] = "PrevMed_Females"
  names(all_res)[names(all_res)=="EstimatePrevF"] = "PrevEst_Females"
  names(all_res)[names(all_res)=="PrevLB_2.5%F"] = "PrevLB_Females"
  names(all_res)[names(all_res)=="PrevUB_97.5%F"] = "PrevUB_Females"

  names(all_res)[names(all_res)=="MedianPrevM"] = "PrevMed_Males"
  names(all_res)[names(all_res)=="EstimatePrevM"] = "PrevEst_Males"
  names(all_res)[names(all_res)=="PrevLB_2.5%M"] = "PrevLB_Males"
  names(all_res)[names(all_res)=="PrevUB_97.5%M"] = "PrevUB_Males"

  names(all_res)[names(all_res)=="MedianPrevM+F"] = "PrevMed_BothSexes"
  names(all_res)[names(all_res)=="EstimatePrevM+F"] = "PrevEst_BothSexes"
  names(all_res)[names(all_res)=="PrevLB_2.5%M+F"] = "PrevLB_BothSexes"
  names(all_res)[names(all_res)=="PrevUB_97.5%M+F"] = "PrevUB_BothSexes"

  #Case prevalence
  names(all_res)[names(all_res)=="CasePrev_EstF"] = "CasePrevEst_Females"
  names(all_res)[names(all_res)=="CasePrev_MedF"] = "CasePrevMed_Females"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%F"] = "CasePrevLB_Females"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%F"] = "CasePrevUB_Females"

  names(all_res)[names(all_res)=="CasePrev_EstM"] = "CasePrevEst_Males"
  names(all_res)[names(all_res)=="CasePrev_MedM"] = "CasePrevMed_Males"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%M"] = "CasePrevLB_Males"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%M"] = "CasePrevUB_Males"

  names(all_res)[names(all_res)=="CasePrev_EstM+F"] = "CasePrevEst_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_MedM+F"] = "CasePrevMed_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_LB_2.5%M+F"] = "CasePrevLB_BothSexes"
  names(all_res)[names(all_res)=="CasePrev_UB_97.5%M+F"] = "CasePrevUB_BothSexes"

  #Incidence rates
  names(all_res)[names(all_res)=="EstimateIncF"] = "InciMed_Females"
  names(all_res)[names(all_res)=="MedianIncF"] = "InciEst_Females"
  names(all_res)[names(all_res)=="IncLB_2.5%F"] = "InciLB_Females"
  names(all_res)[names(all_res)=="IncUB_97.5%F"] = "InciUB_Females"

  names(all_res)[names(all_res)=="EstimateIncM"] = "InciMed_Males"
  names(all_res)[names(all_res)=="MedianIncM"] = "InciEst_Males"
  names(all_res)[names(all_res)=="IncLB_2.5%M"] = "InciLB_Males"
  names(all_res)[names(all_res)=="IncUB_97.5%M"] = "InciUB_Males"

  names(all_res)[names(all_res)=="EstimateIncM+F"] = "InciMed_BothSexes"
  names(all_res)[names(all_res)=="MedianIncM+F"] = "InciEst_BothSexes"
  names(all_res)[names(all_res)=="IncLB_2.5%M+F"] = "InciLB_BothSexes"
  names(all_res)[names(all_res)=="IncUB_97.5%M+F"] = "InciUB_BothSexes"

  #Incidence cases
  names(all_res)[names(all_res)=="CaseInc_EstF"] = "CaseInciMed_Females"
  names(all_res)[names(all_res)=="CaseInc_MedF"] = "CaseInciEst_Females"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%F"] = "CaseInciLB_Females"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%F"] = "CaseInciUB_Females"

  names(all_res)[names(all_res)=="CaseInc_EstM"] = "CaseInciMed_Males"
  names(all_res)[names(all_res)=="CaseInc_MedM"] = "CaseInciEst_Males"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%M"] = "CaseInciLB_Males"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%M"] = "CaseInciUB_Males"

  names(all_res)[names(all_res)=="CaseInc_EstM+F"] = "CaseInciMed_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_MedM+F"] = "CaseInciEst_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_LB_2.5%M+F"] = "CaseInciLB_BothSexes"
  names(all_res)[names(all_res)=="CaseInc_UB_97.5%M+F"] = "CaseInciUB_BothSexes"

  long_all_res <- data.frame();

  for(ctr in unique(all_res$ISO3))#for(ctr in unique(all_res$Country))
  {
    if(!is.element(ctr,all_res$ISO3)) next#if(!is.element(ctr,all_res$Country)) next
    temp_long_ctr <- subset(all_res,ISO3==ctr)#temp_long_ctr <- subset(all_res,Country==ctr)
    iso3  <- ctr
    ctr <- temp_long_ctr$Country[1]

    long_ctr_men <- data.frame(Country = ctr,
                               sex = "Males",
                               datatype = "Model",
                               weight = 1,
                               Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                               indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                               Median = c(temp_long_ctr$CaseInciMed_Males,temp_long_ctr$InciMed_Males,temp_long_ctr$CasePrevMed_Males,temp_long_ctr$PrevMed_Males),
                               BestFit = c(temp_long_ctr$CaseInciEst_Males,temp_long_ctr$InciEst_Males,temp_long_ctr$CasePrevEst_Males,temp_long_ctr$PrevEst_Males),
                               Lower = c(temp_long_ctr$CaseInciLB_Males,temp_long_ctr$InciLB_Males,temp_long_ctr$CasePrevLB_Males,temp_long_ctr$PrevLB_Males),
                               Upper = c(temp_long_ctr$CaseInciUB_Males,temp_long_ctr$InciUB_Males,temp_long_ctr$CasePrevUB_Males,temp_long_ctr$PrevUB_Males),
                               ISO3 = iso3
    )


    long_ctr_women <- data.frame(Country = ctr,
                                 sex = "Females",
                                 datatype = "Model",
                                 weight = 1,
                                 Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                                 indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                                 Median = c(temp_long_ctr$CaseInciMed_Females,temp_long_ctr$InciMed_Females,temp_long_ctr$CasePrevMed_Females,temp_long_ctr$PrevMed_Females),
                                 BestFit = c(temp_long_ctr$CaseInciEst_Females,temp_long_ctr$InciEst_Females,temp_long_ctr$CasePrevEst_Females,temp_long_ctr$PrevEst_Females),
                                 Lower = c(temp_long_ctr$CaseInciLB_Females,temp_long_ctr$InciLB_Females,temp_long_ctr$CasePrevLB_Females,temp_long_ctr$PrevLB_Females),
                                 Upper = c(temp_long_ctr$CaseInciUB_Females,temp_long_ctr$InciUB_Females,temp_long_ctr$CasePrevUB_Females,temp_long_ctr$PrevUB_Females),
                                 ISO3 = iso3
    )

    long_ctr_both <- data.frame(Country = ctr,
                                sex = "BothSexes",
                                datatype = "Model",
                                weight = 1,
                                Year =c(temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year,temp_long_ctr$Year),
                                indicator = rep(c("CaseIncidence","IncidenceRate","CasePrevalence", "PrevalenceRate"), rep(nrow(temp_long_ctr),4)),
                                Median = c(temp_long_ctr$CaseInciMed_BothSexes,temp_long_ctr$InciMed_BothSexes,temp_long_ctr$CasePrevMed_BothSexes,temp_long_ctr$PrevMed_BothSexes),
                                BestFit = c(temp_long_ctr$CaseInciEst_BothSexes,temp_long_ctr$InciEst_BothSexes,temp_long_ctr$CasePrevEst_BothSexes,temp_long_ctr$PrevEst_BothSexes),
                                Lower = c(temp_long_ctr$CaseInciLB_BothSexes,temp_long_ctr$InciLB_BothSexes,temp_long_ctr$CasePrevLB_BothSexes,temp_long_ctr$PrevLB_BothSexes),
                                Upper = c(temp_long_ctr$CaseInciUB_BothSexes,temp_long_ctr$InciUB_BothSexes,temp_long_ctr$CasePrevUB_BothSexes,temp_long_ctr$PrevUB_BothSexes),
                                ISO3 = iso3
    )

    long_ctr <- rbind(long_ctr_men,long_ctr_women,long_ctr_both)

    min_year <- min(long_ctr$Year)
    #Adding data points
    temp_ctr <- subset(SyphData,ISO3_letters==iso3 & Weight_for_Spectrum_fitting>0)#temp_ctr <- subset(SyphData,Country==ctr & Weight_for_Spectrum_fitting>0)
    ppui <- temp_ctr$Prevalence
    ppui[ppui<=0] = 1/100;

    temp_all_res <- data.frame(Country = ctr,
                               sex = NA,
                               datatype = temp_ctr$Data_type,
                               weight = temp_ctr$Weight_for_Spectrum_fitting,
                               Year=temp_ctr$Year,
                               indicator="PrevalenceRate",
                               Median = temp_ctr$Prevalence/100,
                               BestFit = temp_ctr$Prevalence/100,
                               Lower = NA,
                               Upper = NA,
                               ISO3 = iso3
    )

    sdall <- sqrt(ppui/100*(1-ppui/100)/temp_ctr$N_tested)
    sdall[is.na(sdall)] <- 0
    temp_all_res$Lower <- pmax(temp_all_res$Median-1.96*sdall,0)
    temp_all_res$Upper <- pmin(temp_all_res$Median+1.96*sdall,1)

    temp_all_res$sex[temp_all_res$datatype=="ANC Routine screening"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="ANC Survey"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="BloodDonor Screening Men"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="FSW"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="MSM"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="Survey LowRisk Men"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="BloodDonor Screening Men + Women"] <- "BothSexes"
    temp_all_res$sex[temp_all_res$datatype=="Male Sex Workers"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="MSM + MSW combined"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="PWID-Female"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="PWID-Male"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="Survey LowRisk Women"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="Trans-Genders"] <- "BothSexes"
    temp_all_res$sex[temp_all_res$datatype=="BloodDonor Screening Women"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="Prisoners, Men"] <- "Males"
    temp_all_res$sex[temp_all_res$datatype=="Prisoners, Women"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="Wives of PWID"] <- "Females"
    temp_all_res$sex[temp_all_res$datatype=="Survey LowRisk Men+Women"] <- "BothSexes"

    ln_na <- unique(c(which(is.na(temp_all_res$BestFit)),which(is.na(temp_all_res$Median))))
    if(length(ln_na)>0)temp_all_res <- temp_all_res[-ln_na,]
    long_ctr <- rbind(long_ctr,subset(temp_all_res,Year>=min_year))

    long_ctr$datatype <- factor(long_ctr$datatype, levels=c("Model","ANC Routine screening","ANC Survey","BloodDonor Screening Men",
                                                            "FSW","MSM","Survey LowRisk Men","BloodDonor Screening Men + Women",
                                                            "Male Sex Workers","MSM + MSW combined","PWID-Female","PWID-Male",
                                                            "Survey LowRisk Women","Trans-Genders","BloodDonor Screening Women",
                                                            "Prisoners, Men","Prisoners, Women","Wives of PWID","Survey LowRisk Men+Women"))

    #Final merge:
    long_all_res <- rbind(long_all_res,long_ctr)
  }

  m_data <- long_all_res
  plot_aggregate <- function(df, mtitle)
  {
    df$Estimate <- df$BestFit
    df <- subset(df,Year%in%yplot)
    if(nrow(df)==0) return(NULL)
    in_mtitle <- df$Country[1]

    df$datatype <- factor(df$datatype, levels=c("Model","ANC Routine screening","ANC Survey","BloodDonor Screening Men",
                                                "FSW","MSM","Survey LowRisk Men","BloodDonor Screening Men + Women",
                                                "Male Sex Workers","MSM + MSW combined","PWID-Female","PWID-Male",
                                                "Survey LowRisk Women","Trans-Genders","BloodDonor Screening Women",
                                                "Prisoners, Men","Prisoners, Women","Wives of PWID","Survey LowRisk Men+Women"))


    p <- df %>% filter(indicator=="PrevalenceRate")%>% ggplot(subset=.(ISO3==mtitle & datatype=="Model"), aes(x=Year, y=100*Estimate, ymin = 100*Lower, ymax= 100*Upper,color=datatype,fill=datatype))+
        geom_line(data =. %>% filter(ISO3==mtitle& datatype=="Model"), aes(Year, 100*Estimate)) +
        geom_ribbon(data =. %>% filter(ISO3==mtitle& datatype=="Model"),alpha=0.2,linetype=0)+
        geom_pointrange(data =. %>% filter(ISO3==mtitle& datatype!="Model"), position=position_jitter(w = 0.05, h = 0),#position=position_jitter(width=0.5),
                        linetype='solid', aes(size=weight),size=0.75)+
        facet_wrap(~sex)+
        expand_limits(y = 0) +
        labs(title = paste(in_mtitle," Syphilis prevalence trend among adults (15-49 y)", sep=","), x = "Year", y="Test-adjusted prevalence, active syphilis (%)") +
        theme_minimal() +
        theme(legend.position = "bottom")
    p$labels$colour="Source"
    #p$labels$fill="Source"
    return(p)
  }

  by_country <- m_data %>%
    count(Country, BestFit, Median,Lower, Year,Upper, sex, datatype, indicator, ISO3, wt = BestFit, name = "Estimate") %>%
    split(.$ISO3)

  idx <- 1:length(by_country)
  if(!is.null(countries_iso))
  {
    temp_iso <- sapply(1:length(by_country), function(ii){
      res <- "Missing"
      if(!is.null(by_country[[ii]]))
      {
        if(nrow(by_country[[ii]])>=1)
        {
          res <- by_country[[ii]]$ISO3[1]
        }
      }
      return(res)
    })
    idx <- sapply(idx, function(ii){
      if(temp_iso[ii]=="Missing")
        {
          return(-1)
        }else
        {
          if(temp_iso[ii]%in%countries_iso)
          {
            return(ii)
          } else
          {
            return(-1)
          }
        }
      })
  }

  idx <- idx[idx>0]
  allplots <- lapply(idx, function(ii){plot_aggregate(by_country[[ii]],names(by_country)[ii])})
  if(length(idx)==1)
  {
    allplots <- allplots[[1]]
  }
  return(allplots)
}



