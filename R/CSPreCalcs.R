#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
get_precalcsCS <- function(namesinputfiles)
{
  require(openxlsx)

  ################################################################################
  ################################################################################

  fin_name.data.file <- namesinputfiles$prevalence
  fin_name_cs_screening <- namesinputfiles$screening
  fin_name_cs_db <- namesinputfiles$csdb
  fin_name_allbirths <- namesinputfiles$allbirths

  namesCol = c("Country","ISO3" ,"ISO3_letters", "WHO_region","Data_type", "Population_code", "Sex", "Year", "DX_Code",
               "N_positive", "N_tested", "Prevalence", "Weight", "WghtNSpectrum", "Weight_for_Spectrum_fitting")
  TFnamesCol = c("Country","ISO3", "ISO3_letters", "WHO region","Data type", "Population code", "Sex", "Year", "DX_Code",
                 "N positive", "N tested", "Prevalence", "Weight","WghtNSpectrum", "Weight_for_Spectrum_fitting","DiagTestAdjusteFactor")

  high_risk_adj = 1
  zero_prev_adj = 1/100
  require(xlsx)

  SyphData <- openxlsx::read.xlsx(fin_name.data.file,sheet="Data Entry", startRow=1,
                                  cols= 1:(2 + length(namesCol)), check.names = FALSE, sep.names =" ")

  names(SyphData)[names(SyphData)=="Weight"] <- "Weight_for_Spectrum_fitting"
  names(SyphData)[names(SyphData)=="Midpoint study year"] <- "Year"
  names(SyphData)[names(SyphData)=="Population type"] <- "Data type"
  names(SyphData)[names(SyphData)=="Population code"] <- "Data type, code"
  names(SyphData)[names(SyphData)=="DX_Code - updated"] <- "DX_Code"
  names(SyphData)[names(SyphData)=="Duration estimate"] <- "Duration"

  SyphData$ISO3 <- SyphData$ISO3_letters

  SyphData$"Data type" <- sapply(SyphData$"Data type, code", function(xx){
    res <- "Other"
    if(xx==1) res = "ANC Survey" else if(xx==2) res = "ANC Routine screening"
    res
  })


  charnumtransform <- function(x) as.numeric(gsub(",", "", as.character(x)))

  suppressWarnings(SyphData$Prevalence <- charnumtransform(SyphData$Prevalence)) #as.numeric(as.character(SyphData$Prevalence)))
  suppressWarnings(SyphData$`N positive` <- charnumtransform(SyphData$`N positive`))#as.numeric(as.character(SyphData$`N positive`)))
  suppressWarnings(SyphData$`N tested` <- charnumtransform(SyphData$`N tested`))##as.numeric(as.character(SyphData$`N tested`)))
  suppressWarnings(SyphData$Year <- charnumtransform(SyphData$Year))#as.numeric(as.character(SyphData$Year)))

  SyphData <- subset(SyphData, !is.na(Prevalence))
  SyphData$`N tested`[is.na(SyphData$`N tested`) & SyphData$Prevalence==0] <- 300
  SyphData$`N tested`[is.na(SyphData$`N tested`) & is.na(SyphData$`N positive`)] <- 300
  idx_zerprev <- which(SyphData$Prevalence==0)
  SyphData$Prevalence[idx_zerprev] = zero_prev_adj/SyphData$`N tested`[idx_zerprev]

  idx_np <- which(!is.na(SyphData$`N positive`) & (is.na(SyphData$`N tested`)))
  SyphData$`N tested`[idx_np] <- SyphData$`N positive`[idx_np]/SyphData$Prevalence[idx_np] * 100

  idx_nt <- which(is.na(SyphData$`N positive`) & (!is.na(SyphData$`N tested`)))
  SyphData$`N positive`[idx_nt] <- SyphData$`N tested`[idx_nt] * SyphData$Prevalence[idx_nt]/100

  idx_ntp <- which(is.na(SyphData$`N positive`) & (is.na(SyphData$`N tested`)))

  SyphData$`N positive`[idx_ntp] <- 300 * SyphData$Prevalence[idx_ntp]/100
  SyphData$`N tested`[idx_ntp] <- 300

  SyphData$DX_Code[is.na(SyphData$DX_Code)] <- 1
  SyphData <- subset(SyphData, !is.na(SyphData$`N tested`))

  DiagnosticTest <- openxlsx::read.xlsx(fin_name.data.file,sheet="Diagnostic tests", rows=1:8,
                                        cols= 1:6, check.names = FALSE, sep.names =" ")

  colnames(DiagnosticTest)=c("STI name","Diagnostic test", "DX_code", "Adjustment_factor",
                             "Source_for_adjustment_factor", "Comments")

  LowRisk <- c("ANC Routine screening","ANC Survey", "BloodDonor Screening Men", "Survey LowRisk Men", "BloodDonor Screening Men + Women",
               "Survey LowRisk Men+Women", "Survey LowRisk Women", "BloodDonor Screening Women")

  SyphData$Prevalence = sapply(1:length(SyphData$Prevalence),function(ii){
    res = NA;
    if(!is.na(SyphData$DX_Code[ii]))
    {
      hr_adj <- ifelse(any(LowRisk==SyphData$'Data type'[ii]),high_risk_adj,1);
      adj = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==SyphData$DX_Code[ii])];
      res = SyphData$Prevalence[ii]*adj*hr_adj
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

  ISO3 <- SyphDatabase_ISO3
  namesCol = c(namesCol,"DiagTestAdjusteFactor")

  SyphData$ISO3 <- sapply(SyphData$ISO3, function(xx) {
    res <- NA
    if(sum(ISO3$Alpha_ISO3==xx & !is.na(xx) & !is.na(ISO3$Alpha_ISO3))>=1)
    {
      res <- ISO3$Numeric_ISO3[ISO3$Alpha_ISO3==xx & !is.na(xx) & !is.na(ISO3$Alpha_ISO3)][1]
    }
    res
  })

  SyphData = SyphData[, is.element(names(SyphData), TFnamesCol)]
  genv = environment()

  ll = sapply(1:ncol(SyphData), function(jj) if (is.element(names(SyphData)[jj],
                                                            TFnamesCol)) {
    iiind = which(TFnamesCol == names(SyphData)[jj])
    names(genv$SyphData)[jj] = namesCol[iiind]
    jj
  })

  SyphData$WghtNSpectrum = SyphData$Weight_for_Spectrum_fitting
  SyphData <- subset(SyphData, Year >= 2010 & !is.na(Prevalence)& (Data_type == "ANC Routine screening" | Data_type == "ANC Survey"))

  ppui <- SyphData$Prevalence
  ppui[ppui <= 0] = 1/100;
  sdall <- sqrt(ppui/100 * (1 - ppui/100)/SyphData$N_tested)
  SyphData$LowerPrevalence <- pmax(SyphData$Prevalence/100 - 1.96 * sdall, 0)
  SyphData$UpperPrevalence <- pmin(SyphData$Prevalence/100 + 1.96*sdall, 1)
  SyphData$BestPrevalence <- SyphData$Prevalence/100


  SyphDataRaw <- SyphData
  SyphDataRaw <- subset(SyphDataRaw, !is.na(Weight_for_Spectrum_fitting) & (Weight_for_Spectrum_fitting>0))
  #The remove countries with one data point
  all_ctr <- unique(SyphDataRaw$Country)
  num_ctr <- sapply(all_ctr, function(ctr) length(which(SyphDataRaw$Country==ctr)))
  all_ctr_kpt <- all_ctr[num_ctr>=2]
  SyphDataRaw <- subset(SyphDataRaw,Country%in%all_ctr_kpt)

  ####
  SyphDataRaw$SDG_Region <- sapply(SyphDataRaw$ISO3, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  ###############################################################################
  ###############################################################################

  #Preparing Jane data
  JaWbTreatment <- openxlsx::read.xlsx(fin_name_cs_screening,"Treatment")
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="Eastern Mediterranean"] = "EMR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="Africa"] = "AFR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="America"] = "AMR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="Europe"] = "EUR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="South-East Asia"] = "SEAR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="Western Pacific"] = "WPR"

  JaWbTreatment$SDG_Region <- sapply(JaWbTreatment$COUNTRY_CODE, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  JaWb_CongenSyphdb_v0 <- data.frame('Rank, 2012 ABO cases'=rep(NA,nrow(JaWbTreatment)),
                                     'Rank, 2016 ABO cases'=rep(NA,nrow(JaWbTreatment)),
                                     'Rank 2012, CS case RATE'=rep(NA,nrow(JaWbTreatment)),
                                     'Rank 2016, CS case RATE'=rep(NA,nrow(JaWbTreatment)), check.names=FALSE
  )

  JaWb_CongenSyphdb_v0$Country = JaWbTreatment$COUNTRY_NAME
  JaWb_CongenSyphdb_v0$'WHO Region' = JaWbTreatment$REGION_WHO
  JaWb_CongenSyphdb_v0$'ISO code' = JaWbTreatment$COUNTRY_CODE
  JaWb_CongenSyphdb_v0$'ISO3nmb' = NA#JaWbTreatment$ISO3nmb
  JaWb_CongenSyphdb_v0$'Year' = JaWbTreatment$Year
  JaWb_CongenSyphdb_v0$'Live Births' = NA
  JaWb_CongenSyphdb_v0$'Still births' = NA
  JaWb_CongenSyphdb_v0$'Stillbirths, Blencowe & Hogan 2016' = NA
  JaWb_CongenSyphdb_v0$'Pregnancies' = NA
  JaWb_CongenSyphdb_v0$'Still/Live births' = NA
  JaWb_CongenSyphdb_v0$'Source of Live and/or Stillbirths' = NA
  JaWb_CongenSyphdb_v0$'Women with >= 1 ANC visit (%)' = NA #In Screening sheet
  JaWb_CongenSyphdb_v0$'Source of ANC1' = NA #
  JaWb_CongenSyphdb_v0$'N tested (1st ANC visit)' = NA #//Screening Sheet
  JaWb_CongenSyphdb_v0$'N, 1st visits' = NA #//Screening Sheet
  JaWb_CongenSyphdb_v0$'Syphilis-tested (1st ANC, %)' = NA #//Screening Sheet
  JaWb_CongenSyphdb_v0$'N tested (any visit ANC)' = NA #//Screening Sheet
  JaWb_CongenSyphdb_v0$'N any ANC visits' = NA #//Screening Sheet/Spectrum?
  JaWb_CongenSyphdb_v0$'Syphilis-tested (any ANC visit,%)' = NA #
  JaWb_CongenSyphdb_v0$'Source of Test coverage' = NA #
  JaWb_CongenSyphdb_v0$'ANC women with syphilis, treated, N' = JaWbTreatment$'NUM_TREATMENT-TOTAL' #
  JaWb_CongenSyphdb_v0$'Syphilis-infected ANC' = NA #
  JaWb_CongenSyphdb_v0$'Treated (%)' = JaWbTreatment$'PER_TREATMENT-TOTAL' #
  JaWb_CongenSyphdb_v0$'Treated (%)'[JaWb_CongenSyphdb_v0$'Treated (%)'=="A"] <- NA
  JaWb_CongenSyphdb_v0$'Source of Treated' = NA #
  JaWb_CongenSyphdb_v0$'Congenital syphilis case REPORTS' = NA #
  JaWb_CongenSyphdb_v0$'CS case report rate' = NA #

  JaWb_CongenSyphdb_v0$SDG_Region <- sapply(JaWb_CongenSyphdb_v0$'ISO code', function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  #Screening
  ##
  JaWbSCreening <- openxlsx::read.xlsx(fin_name_cs_screening,"Screening", detectDates = TRUE)
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="Eastern Mediterranean"] = "EMR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="Africa"] = "AFR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="America"] = "AMR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="Europe"] = "EUR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="South-East Asia"] = "SEAR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="Western Pacific"] = "WPR"
  names(JaWbSCreening)[names(JaWbSCreening)=="X5"] <- "Coverage"

  JaWbSCreening$SDG_Region <- sapply(JaWbSCreening$COUNTRY_CODE, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  for(ii in 1:nrow(JaWb_CongenSyphdb_v0))
  {
    iso <- JaWb_CongenSyphdb_v0$`ISO code`[ii];
    year_ii <-JaWb_CongenSyphdb_v0$Year[ii]
    tempdata <- subset(JaWbSCreening,COUNTRY_CODE==iso & Year==year_ii)
    if(nrow(tempdata)>=1)
    {
      per_visit <- tempdata$'PER_ANY_VISIT-TOTAL_Calculated'
      per_visit <- per_visit[per_visit!="A" & !is.na(per_visit)]
      if(length(per_visit)>=1)
      {
        per_visit <- mean(as.numeric(as.character(per_visit)))
      } else per_visit <- NA
      JaWb_CongenSyphdb_v0$'Women with >= 1 ANC visit (%)'[ii] = per_visit

      #Number first visits
      num_1visit <- tempdata$'DEN_FIRST_VISIT-TOTAL'
      num_1visit <- num_1visit[num_1visit!="A" & !is.na(num_1visit)]
      if(length(num_1visit)>=1)
      {
        num_1visit <- mean(as.numeric(as.character(num_1visit)))
      } else num_1visit <- NA
      JaWb_CongenSyphdb_v0$'N, 1st visits'[ii] = num_1visit

      #Number tested
      num_tested <- tempdata$'NUM_FIRST_VISIT-TOTAL'
      num_tested <- num_tested[num_tested!="A" &  !is.na(num_tested)]
      if(length(num_tested)>=1)
      {
        num_tested <- mean(as.numeric(as.character(num_tested)))
      } else num_tested <- NA
      JaWb_CongenSyphdb_v0$'N tested (1st ANC visit)'[ii] = num_tested

      #Percentage tested
      per_tested <- tempdata$'PER_FIRST_VISIT-TOTAL_Calculated'
      per_tested <- per_tested[per_tested!="A" &  !is.na(per_tested)]
      if(length(per_tested)>=1)
      {
        per_tested <- mean(as.numeric(as.character(per_tested)))
      } else per_tested <- NA
      JaWb_CongenSyphdb_v0$'Syphilis-tested (1st ANC, %)'[ii] = NA #//Screening Sheet

      #Num tested, any visit
      num_tested <- tempdata$'NUM_ANY_VISIT-TOTAL'
      num_tested <- num_tested[num_tested!="A" &  !is.na(num_tested)]
      if(length(num_tested)>=1)
      {
        num_tested <- mean(as.numeric(as.character(num_tested)))
      } else num_tested <- NA
      JaWb_CongenSyphdb_v0$'N tested (any visit ANC)'[ii] = num_tested

      #Num any visits
      num_tested <- tempdata$'DEN_ANY_VISIT-TOTAL'
      num_tested <- num_tested[num_tested!="A" &  !is.na(num_tested)]
      if(length(num_tested)>=1)
      {
        num_tested <- mean(as.numeric(as.character(num_tested)))
      } else num_tested <- NA
      JaWb_CongenSyphdb_v0$'N any ANC visits'[ii] = num_tested

      #Per any visits
      per_tested <- tempdata$'PER_ANY_VISIT-TOTAL_Calculated'
      per_tested <- per_tested[per_tested!="A" &  !is.na(per_tested)]
      if(length(per_tested)>=1)
      {
        per_tested <- mean(as.numeric(as.character(per_tested)))
      } else per_tested <- NA
      JaWb_CongenSyphdb_v0$'Syphilis-tested (any ANC visit,%)'[ii] = per_tested
    }#End if(nrow(tempdata)>=1)
  }


  JaWb_EarlyANC <- data.frame(ISO3 <- JaWbSCreening$COUNTRY_CODE)
  JaWb_EarlyANC$Country <- JaWbSCreening$COUNTRY_NAME
  JaWb_EarlyANC$Region <- JaWbSCreening$REGION_WHO
  JaWb_EarlyANC$'Start year' <- JaWbSCreening$period_start_date
  JaWb_EarlyANC$'End year' <- JaWbSCreening$period_end_date
  JaWb_EarlyANC$'Coverage early ANC' <- JaWbSCreening$Coverage
  JaWb_EarlyANC$'N early ANC' <- JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Code, first ANC threshold time' <- NA#JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Sample size' <- NA#JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Source code' <- NA#JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Source for Source code' <- NA#JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Early ANC' <- JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Late ANC' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'All women' <- 1#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Before cut-off' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'After cut-off' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'ABO risk, average, treated women' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Time first ANC, before cut-off' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Time first ANC, after cut-off' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Time of first ANC, national average' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Source' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100


  JaWb_EarlyANC$SDG_Region <- sapply(JaWb_EarlyANC$ISO3, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  #########################################################################################################
  #########################################################################################################
  #Reading Eline's data
  ElWb_CongenSyphdb <- openxlsx::read.xlsx(fin_name_cs_db, sheet="CongenSyph db")[,1:30]
  ElWb_CongenSyphdb$COUNTRY[ElWb_CongenSyphdb$COUNTRY=="South-Sudan"] <- "South Sudan"

  names(ElWb_CongenSyphdb) <- c("Rank, 2012 ABO cases", "Rank, 2016 ABO cases", "Rank 2012, CS case RATE", "Rank 2016, CS case RATE", "Country","WHO Region",
                                "ISO code", "ISO3nmb", "Year", "Live Births", "Still births", "Stillbirths, Blencowe & Hogan 2016", "Pregnancies", "Still/Live births",
                                "Source of Live and/or Stillbirths", "Women with >= 1 ANC visit (%)", "Source of ANC1", "N tested (1st ANC visit)", "N, 1st visits",
                                "Syphilis-tested (1st ANC, %)", "N tested (any visit ANC)", "N any ANC visits", "Syphilis-tested (any ANC visit,%)", "Source of Test coverage",
                                "ANC women with syphilis, treated, N", "Syphilis-infected ANC", "Treated (%)", "Source of Treated", "Congenital syphilis case REPORTS",
                                "CS case report rate")

  ElWb_CongenSyphdb$SDG_Region <- sapply(ElWb_CongenSyphdb$'ISO code', function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })


  #Fixing Years
  ElWb_CongenSyphdb <- subset(ElWb_CongenSyphdb,!is.na(Year) & !is.na(ElWb_CongenSyphdb$'ISO code'))
  ElWb_CongenSyphdb$Year[ElWb_CongenSyphdb$Year=="2007-9"] <- "2008"

  ElWb_CongenSyphdb$Year <- as.numeric(as.character(ElWb_CongenSyphdb$Year))
  ElWb_CongenSyphdb$Year[ElWb_CongenSyphdb$Year>3000] <- ElWb_CongenSyphdb$Year[ElWb_CongenSyphdb$Year>=3000]-1000

  #Merging with Jane's data
  ElWb_CongenSyphdb <- rbind(ElWb_CongenSyphdb,JaWb_CongenSyphdb_v0)
  ElWb_CongenSyphdb$imputed_pregnancies <- NA
  ElWb_CongenSyphdb$imputed_ratio_births <- NA
  ElWb_CongenSyphdb$imputed_prevalence <- NA

  Copy_of_ElWb_CongenSyphdb <- ElWb_CongenSyphdb

  ElWb_CongenSyphdb <- ElWb_CongenSyphdb[-c(1:nrow(ElWb_CongenSyphdb)),]

  ######*********************************************************************#####
  #Extend the years
  ######*********************************************************************#####
  ######*
  proj_years <- 1990:2025
  for(tempiso in unique(Copy_of_ElWb_CongenSyphdb$`ISO code`))
  {
    temp_ctr <- subset(Copy_of_ElWb_CongenSyphdb,Copy_of_ElWb_CongenSyphdb$`ISO code`==tempiso)
    if(nrow(temp_ctr)==0) next
    missing_years <- proj_years[which(!is.element(proj_years,temp_ctr$Year))]
    if(length(missing_years)>=1)
    {
      temp_ctr_add <- do.call("rbind", replicate(length(missing_years), temp_ctr[1,], simplify = FALSE))
      temp_ctr_add[,] <- NA
      temp_ctr_add$Year <- missing_years
      temp_ctr_add$Country <- temp_ctr$Country[1]
      temp_ctr_add$`WHO Region` = temp_ctr$`WHO Region`[1]
      temp_ctr_add$`ISO code` <- temp_ctr$`ISO code`[1]
      temp_ctr_add$ISO3nmb <- temp_ctr$ISO3nmb[1]
      temp_ctr <- rbind(temp_ctr,temp_ctr_add)
      rm(temp_ctr_add)
    }

    temp_ctr <- subset(temp_ctr,!is.na(Year))
    idx_sorted_years <- sort.int(temp_ctr$Year,index.return = T)$ix
    temp_ctr <- temp_ctr[idx_sorted_years,]

    temp_ctr$imputed_pregnancies <- ifelse(is.na(temp_ctr$Pregnancies),TRUE,FALSE)
    temp_ctr$imputed_ratio_births <- ifelse(is.na(temp_ctr$imputed_ratio_births),TRUE,FALSE)

    ElWb_CongenSyphdb <- rbind(ElWb_CongenSyphdb,temp_ctr)
  }

  ElWb_CongenSyphdb$'Still births' <- sapply(ElWb_CongenSyphdb$'Still births', function(x)
  {
    res <- NA;
    if(!is.na(x)) res <- as.numeric(as.character(x))
    res
  })

  ElWb_CongenSyphdb$'Live Births' <- sapply(ElWb_CongenSyphdb$'Live Births', function(x)
  {
    res <- NA;
    if(!is.na(x)) res <- as.numeric(as.character(x))
    res
  })

  ElWb_CongenSyphdb$'Still/Live births'[] <- NA;
  ElWb_CongenSyphdb$'Still/Live births' <- ElWb_CongenSyphdb$'Still births'/ElWb_CongenSyphdb$'Live Births'

  #Imputation of missing
  if(any(is.na(ElWb_CongenSyphdb$'Still/Live births')))
  {
    copy_of_ElWb_CongenSyphdb <- ElWb_CongenSyphdb;
    for(ctr in unique(ElWb_CongenSyphdb$`ISO code`))
    {
      idx_ctr <- which(ElWb_CongenSyphdb$`ISO code`==ctr)
      if(length(idx_ctr)==0) next;
      count = 0;
      while((any(is.na(ElWb_CongenSyphdb$'Still/Live births'[idx_ctr]))))
      {
        if(all(is.na(ElWb_CongenSyphdb$'Still/Live births'[idx_ctr]))) #If all are NA,
        {
          #replace with region average
          reg <- ElWb_CongenSyphdb$SDG_Region[idx_ctr[1]]
          #dat_reg <- subset(copy_of_ElWb_CongenSyphdb,copy_of_ElWb_CongenSyphdb$`WHO Region`==reg & !(Country%in%c("AMR total","EMR total","AFR total","SEAR total","WPR total","WPR total")))
          dat_reg <- subset(copy_of_ElWb_CongenSyphdb,copy_of_ElWb_CongenSyphdb$SDG_Region==reg & !(Country%in%c("AMR total","EMR total","AFR total","SEAR total","WPR total","WPR total")))
          ratio <- sapply(ElWb_CongenSyphdb$Year[idx_ctr], function(yy){
            res <- NA;
            rr_yy <- dat_reg$'Still/Live births'[dat_reg$Year==yy]
            if(sum(!is.na(rr_yy))>=1) res <- mean(rr_yy,na.rm=T)
            res
          })
          ElWb_CongenSyphdb$'Still/Live births'[idx_ctr] <- ratio;
        } else if(any(is.na(ElWb_CongenSyphdb$'Still/Live births'[idx_ctr]))) #If only a few are NAs,
        {
          #Interpolate and flatline the edges
          years <- copy_of_ElWb_CongenSyphdb$Year[idx_ctr]
          values <- ElWb_CongenSyphdb$'Still/Live births'[idx_ctr];
          idx_min_na <- which(!is.na(values))
          if(length(idx_min_na)==1)
          {
            ElWb_CongenSyphdb$'Still/Live births'[idx_ctr] = values[idx_min_na];
          } else
          {
            idxsorted <- sort.int(years, index.return=TRUE)$ix
            years <- years[idxsorted];
            values <- values[idxsorted]
            idxmin = min(which(!is.na(values)))[1]
            idxmax = max(which(!is.na(values)))[1]
            if(idxmin > 1) values[1] = values[idxmin]
            if(idxmax<length(values)) values[length(values)] = values[idxmax]
            if(length(values)>2)
            {
              for(i in 2:(length(values)-1))
              {
                if(is.na(values[i]))
                {
                  maxidx <- which(!is.na(values))
                  maxidx <- max(maxidx[maxidx>i])
                  val = values[i-1] + (values[maxidx]-values[i-1])*(years[i]-years[i-1])/(years[maxidx]-years[i-1])
                  values[i] = val

                }
              }
            }#End if(length(values)>2)

            for(i in idx_ctr)
            {
              ElWb_CongenSyphdb$'Still/Live births'[i] = mean(values[which(years==copy_of_ElWb_CongenSyphdb$Year[i])])
            }#End for(i in idx_ctr)
          }
        }#End if(all(is.na(ElWb_CongenSyphdb$'Still/Live births'[idx_ctr]))) #If all are NA
        count = count+1;
        if(count>=2) break
      }#while((any(is.na(ElWb_CongenSyphdb$'Still/Live births'[idx_ctr]))))

    }#End for(ctr in unique(ElWb_CongenSyphdb$Country))
    rm(copy_of_ElWb_CongenSyphdb)
  }

  #options(warn=2)
  #
  ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)' <- sapply(ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)', function(x)
  {
    res <- NA;
    if(!is.na(x))
    {
      cpx <- x
      while(grepl("or",cpx,TRUE))
      {
        idx <- unlist(gregexpr("or",cpx))[1]
        cpx <- substr(cpx,idx+2,nchar(cpx))
      }

      while(grepl("\\?",cpx,TRUE))
      {
        idx <- unlist(gregexpr("\\?",cpx))[1]
        if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
      }
      #cat(x, ", cpx=",cpx,"\n")
      while(grepl("-",cpx,TRUE))
      {
        idx <- unlist(gregexpr("-",cpx))[1]
        if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
      }
      while(grepl("%",cpx,TRUE))
      {
        idx <- unlist(gregexpr("%",cpx))[1]
        if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
      }
      res <- as.numeric(as.character(cpx))
    }
    res
  })

  gm_convertchartonum <- function(cpx)
  {
    while(grepl("or",cpx,TRUE))
    {
      idx <- unlist(gregexpr("or",cpx))[1]
      cpx <- substr(cpx,idx+2,nchar(cpx))
    }

    while(grepl("\\?",cpx,TRUE))
    {
      idx <- unlist(gregexpr("\\?",cpx))[1]
      if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
    }
    #cat(x, ", cpx=",cpx,"\n")
    while(grepl("-",cpx,TRUE))
    {
      idx <- unlist(gregexpr("-",cpx))[1]
      if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
    }
    while(grepl("%",cpx,TRUE))
    {
      idx <- unlist(gregexpr("%",cpx))[1]
      if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
    }

    for(tch in c(">","<",">=", "<=", "="))

      while(grepl(">",cpx,TRUE))
      {
        idx <- unlist(gregexpr(">",cpx))[1]
        if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
      }
    res <- as.numeric(as.character(cpx))
  }

  ####################################
  ###Syphilis-tested (1st ANC, %)
  ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)' <- sapply(ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)', function(x)
  {
    res <- NA;
    if(!is.na(x))
    {
      res <- gm_convertchartonum(x)#as.numeric(as.character(x))
    }
    res
  })

  #Imputation of missing for Syphilis-tested (1st ANC, %)
  if(any(is.na(ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)')))
  {
    copy_of_ElWb_CongenSyphdb <- ElWb_CongenSyphdb;
    for(ctr in sort(unique(ElWb_CongenSyphdb$`ISO code`)))
    {
      idx_ctr <- which(ElWb_CongenSyphdb$`ISO code`==ctr)
      if(length(idx_ctr)==0) next;
      count=0;
      while((any(is.na(ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[idx_ctr]))))
      {
        if(all(is.na(ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[idx_ctr]))) #If all are NA,
        {
          #replace with region average
          #reg <- ElWb_CongenSyphdb$`WHO Region`[idx_ctr[1]]
          dat_reg <- subset(copy_of_ElWb_CongenSyphdb,ElWb_CongenSyphdb$`WHO Region`==reg)
          reg <- ElWb_CongenSyphdb$SDG_Region[idx_ctr[1]]
          dat_reg <- subset(copy_of_ElWb_CongenSyphdb,ElWb_CongenSyphdb$SDG_Region==reg)
          ratio <- sapply(ElWb_CongenSyphdb$Year[idx_ctr], function(yy){
            res <- NA;
            rr_yy <- dat_reg$'Syphilis-tested (1st ANC, %)'[dat_reg$Year==yy]
            if(sum(!is.na(rr_yy))>=1) res <- mean(rr_yy,na.rm=T)
            res
          })
          ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[idx_ctr] <- ratio;
        } else if(any(is.na(ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[idx_ctr]))) #If only a few are NAs,
        {
          #Interpolate and flatline the edges
          years <- copy_of_ElWb_CongenSyphdb$Year[idx_ctr]
          values <- ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[idx_ctr];
          idx_min_na <- which(!is.na(values))

          if(length(idx_min_na)==1)
          {
            ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[idx_ctr] = values[idx_min_na];
          } else
          {
            idxsorted <- sort.int(years, index.return=TRUE)$ix
            years <- years[idxsorted];
            values <- values[idxsorted]
            idxmin = min(which(!is.na(values)))[1]
            idxmax = max(which(!is.na(values)))[1]
            if(idxmin > 1) values[1] = values[idxmin]
            if(idxmax<length(values)) values[length(values)] = values[idxmax]
            if(length(values)>2)
            {
              for(i in 2:(length(values)-1))
              {
                temp_idx <- which(years==years[i])
                if(length(temp_idx)>=2)
                {
                  idx_na <- which(years==years[i] & is.na(values))
                  idx_notna <- which(years==years[i] & !is.na(values))
                  if(length(idx_notna)>=1 & length(idx_na)>=1)
                  {
                    values[idx_na] = mean(values[idx_notna])
                  }
                }#End if(length(temp_idx)>=2)

                if(is.na(values[i]))
                {
                  maxidx <- which(!is.na(values))
                  maxidx <- max(maxidx[maxidx>i])
                  val = values[i-1] + (values[maxidx]-values[i-1])*(years[i]-years[i-1])/(years[maxidx]-years[i-1])
                  values[i] = val
                }#End if(is.na(values[i]))
              }#End for(i in 2:(length(values)-1))
            }#End if(length(values)>2)

            for(i in idx_ctr)
            {
              ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[i] = mean(values[which(years==copy_of_ElWb_CongenSyphdb$Year[i])])
              #cat("i=",i,"\n")
            }#End for(i in idx_ctr)
          }
        }#End if(all(is.na(ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[idx_ctr]))) #If all are NA
        count = count+1
        if(count>=2) break
      }#End while((any(is.na(ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[idx_ctr]))))
    }#End for(ctr in unique(ElWb_CongenSyphdb$Country))
  }

  #####################################
  ###`Treated (%)`
  ElWb_CongenSyphdb$'Treated (%)' <- sapply(ElWb_CongenSyphdb$'Treated (%)', function(x)
  {
    res <- NA;
    if(!is.na(x)& (x!="No positive cases to be treated."))
    {
      res <- gm_convertchartonum(x)#as.numeric(as.character(x))
    }
    res
  })

  #Variable to be converted to numeric
  listvar <- c("Rank, 2012 ABO cases", "Rank, 2016 ABO cases", "Rank 2012, CS case RATE", "Rank 2016, CS case RATE",  #"Country", "WHO Region", "ISO code", "ISO3nmb",
               "Year", "Live Births", "Still births", "Stillbirths, Blencowe & Hogan 2016",
               "Pregnancies", "Still/Live births",  "Women with >= 1 ANC visit (%)", #"Source of ANC1", "Source of Live and/or Stillbirths",
               "N tested (1st ANC visit)", "N, 1st visits", "Syphilis-tested (1st ANC, %)",  "N tested (any visit ANC)", "N any ANC visits",
               "Syphilis-tested (any ANC visit,%)",  "ANC women with syphilis, treated, N", "Syphilis-infected ANC", #"Source of Test coverage",
               "Treated (%)",  "Congenital syphilis case REPORTS", "CS case report rate")#, "Source of Treated",

  for(vname in listvar)
  {
    jj <- which(names(ElWb_CongenSyphdb)==vname)
    ElWb_CongenSyphdb[,jj] <- sapply(ElWb_CongenSyphdb[,jj], function(x)
    {
      res <- NA;
      if(!is.na(x)& (x!="No positive cases to be treated.") & (!grepl("11 of 12 clinics",x)) & (!grepl("Near-universal, near-100",x))&
         (!grepl("https://www.rivm.nl",x))&(!grepl("Numerator is sentinel survey",x))&(!grepl("83%, 64%, 77%, 69%, 76%" ,x))&
         (!grepl("Wijesooriya-S et al." ,x))&(!grepl("Wilkinson-D TMIH" ,x)))
      {
        if(x=="FALSE" | x=="false")
        {
          res <- "0"
        } else if(x=="TRUE"| x=="true")
        {
          res <- "0"
        } else
        {
          res <- gm_convertchartonum(x)#as.numeric(as.character(x))
        }
      }
      res
    })
  }

  #Imputation of missing for Treated (%)
  if(any(is.na(ElWb_CongenSyphdb$'Treated (%)')))
  {
    copy_of_CongenDataIn <- ElWb_CongenSyphdb;
    for(ctr in unique(ElWb_CongenSyphdb$`ISO code`))
    {
      idx_ctr <- which(ElWb_CongenSyphdb$`ISO code`==ctr)
      if(length(idx_ctr)==0) next;

      count =0;
      while((any(is.na(ElWb_CongenSyphdb$'Treated (%)'[idx_ctr]))))
      {
        if(all(is.na(ElWb_CongenSyphdb$'Treated (%)'[idx_ctr]))) #If all are NA,
        {
          #replace with region average
          reg <- ElWb_CongenSyphdb$SDG_Region[idx_ctr[1]]
          dat_reg <- subset(copy_of_CongenDataIn,ElWb_CongenSyphdb$SDG_Region==reg)
          ratio <- sapply(ElWb_CongenSyphdb$Year[idx_ctr], function(yy){
            res <- NA;
            rr_yy <- dat_reg$'Treated (%)'[dat_reg$Year==yy]
            if(sum(!is.na(rr_yy))>=1) res <- mean(rr_yy,na.rm=T)
            res
          })
          ElWb_CongenSyphdb$'Treated (%)'[idx_ctr] <- ratio;
        } else if(any(is.na(ElWb_CongenSyphdb$'Treated (%)'[idx_ctr]))) #If only a few are NAs,
        {
          #Interpolate and flatline the edges
          years <- copy_of_CongenDataIn$Year[idx_ctr]
          values <- ElWb_CongenSyphdb$'Treated (%)'[idx_ctr];
          idx_min_na <- which(!is.na(values))
          if(length(idx_min_na)==1)
          {
            ElWb_CongenSyphdb$'Treated (%)'[idx_ctr] = values[idx_min_na];
          } else
          {
            idxsorted <- sort.int(years, index.return=TRUE)$ix
            years <- years[idxsorted];
            values <- values[idxsorted]
            idxmin = min(which(!is.na(values)))[1]
            idxmax = max(which(!is.na(values)))[1]
            if(idxmin > 1) values[1] = values[idxmin]
            if(idxmax<length(values)) values[length(values)] = values[idxmax]
            if(length(values)>2)
            {
              for(i in 2:(length(values)-1))
              {
                temp_idx <- which(years==years[i])
                if(length(temp_idx)>=2)
                {
                  idx_na <- which(years==years[i] & is.na(values))
                  idx_notna <- which(years==years[i] & !is.na(values))
                  if(length(idx_notna)>=1 & length(idx_na)>=1)
                  {
                    values[idx_na] = mean(values[idx_notna])
                  }
                }#End if(length(temp_idx)>=2)

                if(is.na(values[i]))
                {
                  maxidx <- which(!is.na(values))
                  maxidx <- min(maxidx[maxidx>i])
                  val = values[i-1] + (values[maxidx]-values[i-1])*(years[i]-years[i-1])/(years[maxidx]-years[i-1])
                  values[i] = val
                }#End if(is.na(values[i]))
              }#End for(i in 2:(length(values)-1))
            }#End if(length(values)>2)

            for(i in idx_ctr)
            {
              ElWb_CongenSyphdb$'Treated (%)'[i] = mean(values[which(years==copy_of_CongenDataIn$Year[i])])
            }#End for(i in idx_ctr)
          }
        }#End if(all(is.na(ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)'[idx_ctr]))) #If all are NA

        count = count+1;
        if(count>=2) break;
      }#End while((any(is.na(ElWb_CongenSyphdb$'Treated (%)'[idx_ctr]))))
    }#End for(ctr in unique(ElWb_CongenSyphdb$Country))
    rm(copy_of_CongenDataIn)
  }

  #options(warn=2)
  #***********************NEW*****************************************************#
  #Imputation of missing 'Women with >= 1 ANC visit (%)'
  if(any(is.na(ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)')))
  {
    copy_of_CongenDataIn <- ElWb_CongenSyphdb;
    for(ctr in unique(ElWb_CongenSyphdb$`ISO code`))
    {
      idx_ctr <- which(ElWb_CongenSyphdb$`ISO code`==ctr)
      if(length(idx_ctr)==0) next;

      count =0;
      while((any(is.na(ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)'[idx_ctr]))))
      {
        if(all(is.na(ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)'[idx_ctr]))) #If all are NA,
        {
          #replace with region average
          reg <- ElWb_CongenSyphdb$SDG_Region[idx_ctr[1]]
          dat_reg <- subset(copy_of_CongenDataIn,ElWb_CongenSyphdb$SDG_Region==reg)
          ratio <- sapply(ElWb_CongenSyphdb$Year[idx_ctr], function(yy){
            res <- NA;
            rr_yy <- dat_reg$'Women with >= 1 ANC visit (%)'[dat_reg$Year==yy]
            if(sum(!is.na(rr_yy))>=1) res <- mean(rr_yy,na.rm=T)
            res
          })
          ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)'[idx_ctr] <- ratio;
        } else if(any(is.na(ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)'[idx_ctr]))) #If only a few are NAs,
        {
          #Interpolate and flatline the edges
          years <- copy_of_CongenDataIn$Year[idx_ctr]
          values <- ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)'[idx_ctr];
          idx_min_na <- which(!is.na(values))
          if(length(idx_min_na)==1)
          {
            ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)'[idx_ctr] = values[idx_min_na];
          } else
          {
            idxsorted <- sort.int(years, index.return=TRUE)$ix
            years <- years[idxsorted];
            values <- values[idxsorted]
            idxmin = min(which(!is.na(values)))[1]
            idxmax = max(which(!is.na(values)))[1]
            if(idxmin > 1) values[1] = values[idxmin]
            if(idxmax<length(values)) values[length(values)] = values[idxmax]
            if(length(values)>2)
            {
              for(i in 2:(length(values)-1))
              {
                temp_idx <- which(years==years[i])
                if(length(temp_idx)>=2)
                {
                  #kk
                  idx_na <- which(years==years[i] & is.na(values))
                  idx_notna <- which(years==years[i] & !is.na(values))
                  if(length(idx_notna)>=1 & length(idx_na)>=1)
                  {
                    values[idx_na] = mean(values[idx_notna])
                  }
                }#End if(length(temp_idx)>=2)

                if(is.na(values[i]))
                {
                  maxidx <- which(!is.na(values))
                  maxidx <- min(maxidx[maxidx>i])
                  val = values[i-1] + (values[maxidx]-values[i-1])*(years[i]-years[i-1])/(years[maxidx]-years[i-1])
                  values[i] = val
                }#End if(is.na(values[i]))
              }#End for(i in 2:(length(values)-1))
            }#End if(length(values)>2)

            for(i in idx_ctr)
            {
              ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)'[i] = mean(values[which(years==copy_of_CongenDataIn$Year[i])])
            }#End for(i in idx_ctr)
          }
        }#End if(all(is.na(ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)'[idx_ctr]))) #If all are NA

        count = count+1;
        if(count>=2) break;
      }#End while((any(is.na(ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)'[idx_ctr]))))
    }#End for(ctr in unique(ElWb_CongenSyphdb$Country))
  }

  #***********************End NEW*****************************************************#

  ################################################################################
  ######Imputation of pregnancies#################################################
  ################################################################################
  SpectrumBirths <- openxlsx::read.xlsx(fin_name_allbirths, startRow = 2)
  #Countries with missing Spectrum files, we interpolate using a linear model
  for(ctr in unique(ElWb_CongenSyphdb$Country))
  {
    idx_ctr <- which(ElWb_CongenSyphdb$Country==ctr)
    if(any(is.na(ElWb_CongenSyphdb$Pregnancies[idx_ctr])))
    {
      idx_na <-  which(is.na(ElWb_CongenSyphdb$Pregnancies) & ElWb_CongenSyphdb$Country==ctr & !is.na(ElWb_CongenSyphdb$Year))
      idx_notna <-  which(!is.na(ElWb_CongenSyphdb$Pregnancies) & ElWb_CongenSyphdb$Country==ctr & !is.na(ElWb_CongenSyphdb$Year))
      if(length(idx_notna)>=2 & length(idx_na)>=1)
      {
        mod <- lm(Pregnancies~Year,data=ElWb_CongenSyphdb[idx_notna,])
        pred <- predict(mod, newdata=ElWb_CongenSyphdb[idx_na,])
        pred[pred<0] <- 0
        ElWb_CongenSyphdb$Pregnancies[idx_na] <- pred
      }
    }#End if(any(is.na(ElWb_CongenSyphdb$Pregnancies[idx_ctr])))
  }


  for(ii in which(is.na(ElWb_CongenSyphdb$Pregnancies) | is.na(ElWb_CongenSyphdb$`Live Births`)))
  {
    iso = ElWb_CongenSyphdb$`ISO code`[ii]
    year_ii <- ElWb_CongenSyphdb$Year[ii]
    if(is.na(year_ii)) next
    ratiostilllive <- ElWb_CongenSyphdb$`Still/Live births`[ii]
    livebirths <- SpectrumBirths[SpectrumBirths$ISO3==iso, names(SpectrumBirths) == as.character(year_ii)]
    if(length(livebirths)>=1)
    {
      livebirths <- sum(livebirths)
      still <- ratiostilllive/(1+ratiostilllive)*livebirths;
      live <- livebirths/(1+ratiostilllive);
      ElWb_CongenSyphdb$`Live Births`[ii] <- live;
      ElWb_CongenSyphdb$`Still births`[ii] <- still;
      ElWb_CongenSyphdb$Pregnancies[ii] <- still+live;
    } else
    {
      livebirths <- ElWb_CongenSyphdb$Pregnancies[ii]
      still <- ratiostilllive/(1+ratiostilllive)*livebirths;
      live <- livebirths/(1+ratiostilllive);
      ElWb_CongenSyphdb$`Live Births`[ii] <- live;
      ElWb_CongenSyphdb$`Still births`[ii] <- still;
    }
  }

  ElWb_CongenSyphdb$SDG_Region <- sapply(ElWb_CongenSyphdb$'ISO code', function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })
  ################################################################################
  ################################################################################
  ################################################################################
  ElWb_RiskABO_Asumptions <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Assumptions",rows=2:5, cols=1:10)
  names(ElWb_RiskABO_Asumptions)<- c("X1","treated mothers", "CS, risk: liveborn", "CS risk: stillbirth", "CS risk: neonatal death",
                                     "CS risk: prematurity or LBW", "All outcomes", "LB on ABO risk assumption", "UB on ABO risk assumption",
                                     "Variance on ABO risk assumption")

  ElWb_SE_Asumptions <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Assumptions",rows=11:17, cols=1:3)
  ElWb_ABO_Asumptions <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Assumptions",rows=10:17, cols=10:14)
  names(ElWb_ABO_Asumptions) <- c("WHO Regions", "Indicator","ABO risk, treated mothers","Average timing of first ANC", "Number of Pregnancies, 2016")
  #
  ElWb_TimingANC_Asumptions <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Assumptions",rows=20:30, cols=1:10)
  names(ElWb_TimingANC_Asumptions) <- c("Code, in data", "Timing of first ANC", "Median weeks in the pregnancy, of first ANC", "ABO risk probability/treated syphilis-infected mother",
                                        "Minimum (weeks)", "Maximum (weeks)", "Average (weeks)", "ABO risk probability/untreated syphilis-infected mother", "Added risk, from early to late",
                                        "Frequency")

  ElWb_EarlyANC <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Early ANC",rows=1:528, cols=1:21)
  names(ElWb_EarlyANC) <- c("ISO3", "Country", "Region", "Start year", "End year", "Coverage early ANC", "N early ANC", "Code, first ANC threshold time",
                            "Sample size", "Source code", "Source for Source code", "Early ANC", "Late ANC", "All women", "Before cut-off", "After cut-off",
                            "ABO risk, average, treated women", "Time first ANC, before cut-off", "Time first ANC, after cut-off",
                            "Time of first ANC, national average", "Source")


  ElWb_EarlyANC$SDG_Region <- sapply(ElWb_EarlyANC$'ISO3', function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  ElWb_EarlyANC$`Start year`[is.na(ElWb_EarlyANC$`Start year`)] <- ElWb_EarlyANC$`End year`[is.na(ElWb_EarlyANC$`Start year`)]
  ElWb_EarlyANC <- subset(ElWb_EarlyANC, (!is.na(ElWb_EarlyANC$`Start year`)) & (!is.na(ElWb_EarlyANC$`End year`)))


  LongEarlyANCData <- data.frame();
  for(ctriso in unique(ElWb_EarlyANC$ISO3))
  {
    if(is.na(ctriso)) next;
    long_ctr_data <- data.frame()
    ctr_data <- subset(ElWb_EarlyANC,ISO3==ctriso)
    for(ii in 1:nrow(ctr_data))
    {
      ctr_dat_slice <- ctr_data[ii,]
      n_years <- ctr_dat_slice$`End year`[1]-ctr_dat_slice$`Start year`[1]+1
      if(n_years<0)
      {
        n_years = -n_years;
        fy <- ctr_dat_slice$`Start year`[1]
        ctr_dat_slice$`Start year`[1] =ctr_dat_slice$`End year`[1]
        ctr_dat_slice$`End year`[1] = fy
      }
      temp_long_ctr_data <- do.call("rbind", replicate(n_years, ctr_dat_slice, simplify = FALSE))
      temp_long_ctr_data$`Start year` = ctr_dat_slice$`Start year`[1]+0:(n_years-1)
      temp_long_ctr_data$`End year` = ctr_dat_slice$`Start year`[1]+1:n_years;
      long_ctr_data <- rbind(long_ctr_data,temp_long_ctr_data)
    }

    nanames <- c("Coverage early ANC","N early ANC","Code, first ANC threshold time", "Sample size",
                 "Source code","Source for Source code","Early ANC","Late ANC","All women","Before cut-off",
                 "After cut-off","ABO risk, average, treated women","Time first ANC, before cut-off","Time first ANC, after cut-off",
                 "Time of first ANC, national average","Source","SDG_Region")

    temp_proj_years <- 1980:2025
    fin_long_ctr_data <- data.frame();
    for(yy in temp_proj_years)
    {
      temp_long_ctr_data <- subset(long_ctr_data,long_ctr_data$'Start year'==yy)
      if(nrow(temp_long_ctr_data)==0)
      {
        temp_long_ctr_data <- long_ctr_data[1,]
        temp_long_ctr_data$`Start year` = yy
        temp_long_ctr_data$`End year` = yy+1
        temp_long_ctr_data[,is.element(names(temp_long_ctr_data),nanames)] <- NA
      } else if (nrow(temp_long_ctr_data)>1)
      {
        temp_long_ctr_data0 <- temp_long_ctr_data
        temp_long_ctr_data <- temp_long_ctr_data[1,]
        for(vname in nanames[c(1,2,4,7,8,9,10,11,12,13,14,15)])
        {
          temp_long_ctr_data[,names(temp_long_ctr_data)==vname] <- mean(temp_long_ctr_data0[,names(temp_long_ctr_data0)==vname],na.rm=T)
        }
      }
      fin_long_ctr_data <- rbind(fin_long_ctr_data,temp_long_ctr_data)
    }

    #Interpolate and flatline the edges
    for(vname in nanames[c(1,2,4,7,8,9,10,11,12,13,14,15)])
    {
      idx_var <- which(names(fin_long_ctr_data)==vname)
      values <- fin_long_ctr_data[,idx_var]
      if(is.character(values)) next;
      idx_min_na <- which(!is.na(values))
      if(length(idx_min_na)==0) next
      if(length(idx_min_na)==1)
      {
        values <- values[idx_min_na]
      } else if(length(idx_min_na)!=length(values))
      {
        idxmin = min(which(!is.na(values)))[1]
        idxmax = max(which(!is.na(values)))[1]
        if(idxmin > 1) values[1] = values[idxmin]
        if(idxmax<length(values)) values[length(values)] = values[idxmax]
        if(length(values)>2)
        {
          for(i in 2:(length(values)-1))
          {
            if(is.na(values[i]))
            {
              maxidx <- which(!is.na(values))
              maxidx <- min(maxidx[maxidx>i])
              val = values[i-1] + (values[maxidx]-values[i-1])*(fin_long_ctr_data$`Start year`[i]-fin_long_ctr_data$`Start year`[i-1])/(fin_long_ctr_data$`Start year`[maxidx]-fin_long_ctr_data$`Start year`[i-1])
              values[i] = val
            }#End if(is.na(values[i]))
          }#End for(i in 2:(length(values)-1))
        }#End if(length(values)>2)
      }
      fin_long_ctr_data[,idx_var] <- values
    }
    LongEarlyANCData <- rbind(LongEarlyANCData,fin_long_ctr_data)
  }


  ################################################################################
  ################################################################################
  fn_impute <- function(vname,dftoimpute, dfcandidate)
  {
    idx_var <- which(names(dftoimpute)==vname)
    values <- dftoimpute[,idx_var]
    years <- dftoimpute$Year
    if(is.character(values)) return(values);
    if(all(is.na(values)))
    {
      idx_varcandiate <- which(names(dfcandidate)==vname)
      reg <- dftoimpute$SDG_Region[1]
      temp1 <- subset(dfcandidate,(dfcandidate$SDG_Region==reg))
      if(all(is.na(temp1[,idx_varcandiate])))
      {
        temp1 <- dfcandidate;
      }#End if(all(is.na(temp1[,idx_varcandiate])))

      for(ii in seq_len(length(values)))
      {
        #print(ii)
        temp11 <- subset(temp1,Year==years[ii])
        if(nrow(temp11)>=1)
        {
          if((!all(is.na(temp11[,idx_varcandiate])))) values[ii] <- mean(temp11[,idx_varcandiate],na.rm=T)
        }
      }#End for(ii in 1:seq_len(length(values)))
    }#End if(all(is.na(values)))

    if(all(is.na(values))) return(values)

    idxsorted <- sort.int(years, index.return=TRUE)$ix
    sortedyears <- years[idxsorted];
    sortedvalues <- values[idxsorted]
    idxmin = min(which(!is.na(sortedvalues)))[1]
    idxmax = max(which(!is.na(sortedvalues)))[1]
    if(idxmin > 1) sortedvalues[1] = sortedvalues[idxmin]
    if(idxmax<length(sortedvalues)) sortedvalues[length(sortedvalues)] = sortedvalues[idxmax]
    if(length(sortedvalues)>2)
    {
      for(i in 2:(length(sortedvalues)-1))
      {
        temp_idx <- which(sortedyears==years[i])
        if(length(temp_idx)>=2)
        {
          idx_na <- which(sortedyears==sortedyears[i] & is.na(sortedvalues))
          idx_notna <- which(sortedyears==sortedyears[i] & !is.na(sortedvalues))
          if(length(idx_notna)>=1 & length(idx_na)>=1)
          {
            sortedvalues[idx_na] = mean(sortedvalues[idx_notna])
          }
        }#End if(length(temp_idx)>=2)

        if(is.na(sortedvalues[i]))
        {
          maxidx <- which(!is.na(sortedvalues))
          maxidx <- max(maxidx[maxidx>i])
          val = sortedvalues[i-1] + (sortedvalues[maxidx]-sortedvalues[i-1])*(sortedyears[i]-sortedyears[i-1])/(sortedyears[maxidx]-sortedyears[i-1])
          sortedvalues[i] = val
        }#End if(is.na(sortedvalues[i]))
      }#End for(i in 2:(length(sortedvalues)-1))
    }#End if(length(values)>2)

    for(i in 1:length(sortedvalues))
    {
      values[i] = mean(sortedvalues[which(sortedyears==years[i])])
    }#End for(i in 1:length(sortedvalues))
    return(values)
  }

  CongenDataIn <- ElWb_CongenSyphdb; rm(ElWb_CongenSyphdb)
  CongenDataIn$Country[CongenDataIn$Country=="Saint Vincent and Grenadines"] = "Saint Vincent and the Grenadines"
  CongenDataIn$Country[CongenDataIn$Country=="Bahamas (the)"] = "Bahamas"
  CongenDataIn$Country[CongenDataIn$Country=="Cabo Verde"] = "Cape Verde"
  CongenDataIn$Country[CongenDataIn$Country=="Central African Republic (the)"] = "Central African Republic"
  CongenDataIn$Country[CongenDataIn$Country=="Comoros (the)"] = "Comoros"
  cst = CongenDataIn$Country[grep("Brazzaville", CongenDataIn$Country)][1]
  if(length(cst)>=1)
  {
    CongenDataIn$Country[CongenDataIn$Country==cst] = "Congo"
  }
  CongenDataIn$Country[CongenDataIn$Country=="Congo (Brazzaville)"] = "Congo"
  CongenDataIn$Country[CongenDataIn$Country=="C?te d'Ivoire"] = "Cote d'Ivoire"
  CongenDataIn$Country[CongenDataIn$Country=="Dem. People's Republic of Korea"] = "Democratic People's Republic of Korea"
  CongenDataIn$Country[CongenDataIn$Country=="Dominican Republic (the)"] = "Dominican Republic"
  CongenDataIn$Country[CongenDataIn$Country=="Swaziland"] = "Eswatini"
  CongenDataIn$Country[CongenDataIn$Country=="Micronesia (Fed. States of)"] = "Micronesia (Federated States of)"
  CongenDataIn$Country[CongenDataIn$Country=="The Gambia"] = "Gambia"
  CongenDataIn$Country[CongenDataIn$Country=="United Kingdom of Great Britain and Northern Ireland"] = "United Kingdom"
  CongenDataIn$Country[CongenDataIn$Country=="Venezuela, Bolivarian Republic of"] = "Venezuela (Bolivarian Republic of)"
  CongenDataIn$Country[CongenDataIn$Country=="Viet Nam"] = "Vietnam"
  CongenDataIn$Country[CongenDataIn$Country=="Czech Republic"] = "Czechia"

  RiskABO_Asumptions <- ElWb_RiskABO_Asumptions; rm(ElWb_RiskABO_Asumptions)
  SE_Asumptions <- ElWb_SE_Asumptions; rm(ElWb_SE_Asumptions)
  TimingANC_Asumptions <- ElWb_TimingANC_Asumptions; rm(ElWb_TimingANC_Asumptions)
  EarlyANC <- ElWb_EarlyANC; rm(ElWb_EarlyANC)

  result <- list(CongenDataIn=CongenDataIn,SDGRegions=SDGRegions,RiskABO_Asumptions=RiskABO_Asumptions,
                   LongEarlyANCData=LongEarlyANCData, SE_Asumptions=SE_Asumptions,TimingANC_Asumptions=TimingANC_Asumptions,EarlyANC=EarlyANC,
                   fn_impute=fn_impute)
  return(result)
}

#' A function rather aimed at developers
#' @description A function that does blabla, blabla.
#' @keywords internal
#' @export
get_rawCS <- function(namesinputfiles)
{
  require(openxlsx)

  ################################################################################
  ################################################################################

  fin_name.data.file <- namesinputfiles$prevalence
  fin_name_cs_screening <- namesinputfiles$screening
  fin_name_cs_db <- namesinputfiles$csdb
  fin_name_allbirths <- namesinputfiles$allbirths

  namesCol = c("Country","ISO3" ,"ISO3_letters", "WHO_region","Data_type", "Population_code", "Sex", "Year", "DX_Code",
               "N_positive", "N_tested", "Prevalence", "Weight", "WghtNSpectrum", "Weight_for_Spectrum_fitting")
  TFnamesCol = c("Country","ISO3", "ISO3_letters", "WHO region","Data type", "Population code", "Sex", "Year", "DX_Code",
                 "N positive", "N tested", "Prevalence", "Weight","WghtNSpectrum", "Weight_for_Spectrum_fitting","DiagTestAdjusteFactor")

  high_risk_adj = 1
  zero_prev_adj = 1/100
  require(xlsx)

  SyphData <- openxlsx::read.xlsx(fin_name.data.file,sheet="Data Entry", startRow=1,
                                  cols= 1:(2 + length(namesCol)), check.names = FALSE, sep.names =" ")

  names(SyphData)[names(SyphData)=="Weight"] <- "Weight_for_Spectrum_fitting"
  names(SyphData)[names(SyphData)=="Midpoint study year"] <- "Year"
  names(SyphData)[names(SyphData)=="Population type"] <- "Data type"
  names(SyphData)[names(SyphData)=="Population code"] <- "Data type, code"
  names(SyphData)[names(SyphData)=="DX_Code - updated"] <- "DX_Code"
  names(SyphData)[names(SyphData)=="Duration estimate"] <- "Duration"

  SyphData$ISO3 <- SyphData$ISO3_letters

  SyphData$"Data type" <- sapply(SyphData$"Data type, code", function(xx){
    res <- "Other"
    if(xx==1) res = "ANC Survey" else if(xx==2) res = "ANC Routine screening"
    res
  })

  charnumtransform <- function(x) as.numeric(gsub(",", "", as.character(x)))

  suppressWarnings(SyphData$Prevalence <- charnumtransform(SyphData$Prevalence)) #as.numeric(as.character(SyphData$Prevalence)))
  suppressWarnings(SyphData$`N positive` <- charnumtransform(SyphData$`N positive`)) #as.numeric(as.character(SyphData$`N positive`)))
  suppressWarnings(SyphData$`N tested` <- charnumtransform(SyphData$`N tested`)) #as.numeric(as.character(SyphData$`N tested`)))
  suppressWarnings(SyphData$Year <- charnumtransform(SyphData$Year)) #as.numeric(as.character(SyphData$Year)))

  SyphData <- subset(SyphData, !is.na(Prevalence))
  SyphData$`N tested`[is.na(SyphData$`N tested`) & SyphData$Prevalence==0] <- 300
  SyphData$`N tested`[is.na(SyphData$`N tested`) & is.na(SyphData$`N positive`)] <- 300
  idx_zerprev <- which(SyphData$Prevalence==0)
  SyphData$Prevalence[idx_zerprev] = zero_prev_adj/SyphData$`N tested`[idx_zerprev]

  idx_np <- which(!is.na(SyphData$`N positive`) & (is.na(SyphData$`N tested`)))
  SyphData$`N tested`[idx_np] <- SyphData$`N positive`[idx_np]/SyphData$Prevalence[idx_np] * 100

  idx_nt <- which(is.na(SyphData$`N positive`) & (!is.na(SyphData$`N tested`)))
  SyphData$`N positive`[idx_nt] <- SyphData$`N tested`[idx_nt] * SyphData$Prevalence[idx_nt]/100

  idx_ntp <- which(is.na(SyphData$`N positive`) & (is.na(SyphData$`N tested`)))

  SyphData$`N positive`[idx_ntp] <- 300 * SyphData$Prevalence[idx_ntp]/100
  SyphData$`N tested`[idx_ntp] <- 300

  SyphData$DX_Code[is.na(SyphData$DX_Code)] <- 1
  SyphData <- subset(SyphData, !is.na(SyphData$`N tested`))

  DiagnosticTest <- openxlsx::read.xlsx(fin_name.data.file,sheet="Diagnostic tests", rows=1:8,
                                        cols= 1:6, check.names = FALSE, sep.names =" ")

  colnames(DiagnosticTest)=c("STI name","Diagnostic test", "DX_code", "Adjustment_factor",
                             "Source_for_adjustment_factor", "Comments")

  LowRisk <- c("ANC Routine screening","ANC Survey", "BloodDonor Screening Men", "Survey LowRisk Men", "BloodDonor Screening Men + Women",
               "Survey LowRisk Men+Women", "Survey LowRisk Women", "BloodDonor Screening Women")

  SyphData$Prevalence = sapply(1:length(SyphData$Prevalence),function(ii){
    res = NA;
    if(!is.na(SyphData$DX_Code[ii]))
    {
      hr_adj <- ifelse(any(LowRisk==SyphData$'Data type'[ii]),high_risk_adj,1);
      adj = DiagnosticTest$Adjustment_factor[which(DiagnosticTest$DX_code==SyphData$DX_Code[ii])];
      res = SyphData$Prevalence[ii]*adj*hr_adj
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


  ISO3 <- SyphDatabase_ISO3

  namesCol = c(namesCol,"DiagTestAdjusteFactor")

  SyphData$ISO3 <- sapply(SyphData$ISO3, function(xx) {
    res <- NA
    if(sum(ISO3$Alpha_ISO3==xx & !is.na(xx) & !is.na(ISO3$Alpha_ISO3))>=1)
    {
      res <- ISO3$Numeric_ISO3[ISO3$Alpha_ISO3==xx & !is.na(xx) & !is.na(ISO3$Alpha_ISO3)][1]
    }
    res
  })

  SyphData = SyphData[, is.element(names(SyphData), TFnamesCol)]
  genv = environment()

  ll = sapply(1:ncol(SyphData), function(jj) if (is.element(names(SyphData)[jj],
                                                            TFnamesCol)) {
    iiind = which(TFnamesCol == names(SyphData)[jj])
    names(genv$SyphData)[jj] = namesCol[iiind]
    jj
  })

  SyphData$WghtNSpectrum = SyphData$Weight_for_Spectrum_fitting
  SyphData <- subset(SyphData, Year >= 2010 & !is.na(Prevalence)& (Data_type == "ANC Routine screening" | Data_type == "ANC Survey"))

  ppui <- SyphData$Prevalence
  ppui[ppui <= 0] = 1/100;
  sdall <- sqrt(ppui/100 * (1 - ppui/100)/SyphData$N_tested)
  SyphData$LowerPrevalence <- pmax(SyphData$Prevalence/100 - 1.96 * sdall, 0)
  SyphData$UpperPrevalence <- pmin(SyphData$Prevalence/100 + 1.96*sdall, 1)
  SyphData$BestPrevalence <- SyphData$Prevalence/100


  SyphDataRaw <- SyphData
  SyphDataRaw <- subset(SyphDataRaw, !is.na(Weight_for_Spectrum_fitting) & (Weight_for_Spectrum_fitting>0))
  #The remove countries with one data point
  all_ctr <- unique(SyphDataRaw$Country)
  num_ctr <- sapply(all_ctr, function(ctr) length(which(SyphDataRaw$Country==ctr)))
  all_ctr_kpt <- all_ctr[num_ctr>=2]
  SyphDataRaw <- subset(SyphDataRaw,Country%in%all_ctr_kpt)

  ####
  SyphDataRaw$SDG_Region <- sapply(SyphDataRaw$ISO3, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })


  ###############################################################################
  ###############################################################################
  #Preparing Jane's data
  JaWbTreatment <- openxlsx::read.xlsx(fin_name_cs_screening,"Treatment")
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="Eastern Mediterranean"] = "EMR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="Africa"] = "AFR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="America"] = "AMR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="Europe"] = "EUR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="South-East Asia"] = "SEAR"
  JaWbTreatment$REGION_WHO[JaWbTreatment$REGION_WHO=="Western Pacific"] = "WPR"


  JaWbTreatment$SDG_Region <- sapply(JaWbTreatment$COUNTRY_CODE, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  JaWb_CongenSyphdb_v0 <- data.frame('Rank, 2012 ABO cases'=rep(NA,nrow(JaWbTreatment)),
                                     'Rank, 2016 ABO cases'=rep(NA,nrow(JaWbTreatment)),
                                     'Rank 2012, CS case RATE'=rep(NA,nrow(JaWbTreatment)),
                                     'Rank 2016, CS case RATE'=rep(NA,nrow(JaWbTreatment)), check.names=FALSE
  )

  JaWb_CongenSyphdb_v0$Country = JaWbTreatment$COUNTRY_NAME
  JaWb_CongenSyphdb_v0$'WHO Region' = JaWbTreatment$REGION_WHO
  JaWb_CongenSyphdb_v0$'ISO code' = JaWbTreatment$COUNTRY_CODE
  JaWb_CongenSyphdb_v0$'ISO3nmb' = NA#JaWbTreatment$ISO3nmb
  JaWb_CongenSyphdb_v0$'Year' = JaWbTreatment$Year
  JaWb_CongenSyphdb_v0$'Live Births' = NA
  JaWb_CongenSyphdb_v0$'Still births' = NA
  JaWb_CongenSyphdb_v0$'Stillbirths, Blencowe & Hogan 2016' = NA
  JaWb_CongenSyphdb_v0$'Pregnancies' = NA
  JaWb_CongenSyphdb_v0$'Still/Live births' = NA
  JaWb_CongenSyphdb_v0$'Source of Live and/or Stillbirths' = NA
  JaWb_CongenSyphdb_v0$'Women with >= 1 ANC visit (%)' = NA #In Screening sheet
  JaWb_CongenSyphdb_v0$'Source of ANC1' = NA #
  JaWb_CongenSyphdb_v0$'N tested (1st ANC visit)' = NA #//Screening Sheet
  JaWb_CongenSyphdb_v0$'N, 1st visits' = NA #//Screening Sheet
  JaWb_CongenSyphdb_v0$'Syphilis-tested (1st ANC, %)' = NA #//Screening Sheet
  JaWb_CongenSyphdb_v0$'N tested (any visit ANC)' = NA #//Screening Sheet
  JaWb_CongenSyphdb_v0$'N any ANC visits' = NA #//Screening Sheet/Spectrum?
  JaWb_CongenSyphdb_v0$'Syphilis-tested (any ANC visit,%)' = NA #
  JaWb_CongenSyphdb_v0$'Source of Test coverage' = NA #
  JaWb_CongenSyphdb_v0$'ANC women with syphilis, treated, N' = JaWbTreatment$'NUM_TREATMENT-TOTAL' #
  JaWb_CongenSyphdb_v0$'Syphilis-infected ANC' = NA #
  JaWb_CongenSyphdb_v0$'Treated (%)' = JaWbTreatment$'PER_TREATMENT-TOTAL' #
  JaWb_CongenSyphdb_v0$'Treated (%)'[JaWb_CongenSyphdb_v0$'Treated (%)'=="A"] <- NA
  JaWb_CongenSyphdb_v0$'Source of Treated' = NA #
  JaWb_CongenSyphdb_v0$'Congenital syphilis case REPORTS' = NA #
  JaWb_CongenSyphdb_v0$'CS case report rate' = NA #

  JaWb_CongenSyphdb_v0$SDG_Region <- sapply(JaWb_CongenSyphdb_v0$'ISO code', function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  #Screening
  ##
  JaWbSCreening <- openxlsx::read.xlsx(fin_name_cs_screening,"Screening", detectDates = TRUE)
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="Eastern Mediterranean"] = "EMR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="Africa"] = "AFR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="America"] = "AMR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="Europe"] = "EUR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="South-East Asia"] = "SEAR"
  JaWbSCreening$REGION_WHO[JaWbSCreening$REGION_WHO=="Western Pacific"] = "WPR"
  names(JaWbSCreening)[names(JaWbSCreening)=="X5"] <- "Coverage"

  JaWbSCreening$SDG_Region <- sapply(JaWbSCreening$COUNTRY_CODE, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  for(ii in 1:nrow(JaWb_CongenSyphdb_v0))
  {
    iso <- JaWb_CongenSyphdb_v0$`ISO code`[ii];
    year_ii <-JaWb_CongenSyphdb_v0$Year[ii]
    tempdata <- subset(JaWbSCreening,COUNTRY_CODE==iso & Year==year_ii)
    if(nrow(tempdata)>=1)
    {
      per_visit <- tempdata$'PER_ANY_VISIT-TOTAL_Calculated'
      per_visit <- per_visit[per_visit!="A" & !is.na(per_visit)]
      if(length(per_visit)>=1)
      {
        per_visit <- mean(as.numeric(as.character(per_visit)))
      } else per_visit <- NA
      JaWb_CongenSyphdb_v0$'Women with >= 1 ANC visit (%)'[ii] = per_visit

      #Number first visits
      num_1visit <- tempdata$'DEN_FIRST_VISIT-TOTAL'
      num_1visit <- num_1visit[num_1visit!="A" & !is.na(num_1visit)]
      if(length(num_1visit)>=1)
      {
        num_1visit <- mean(as.numeric(as.character(num_1visit)))
      } else num_1visit <- NA
      JaWb_CongenSyphdb_v0$'N, 1st visits'[ii] = num_1visit

      #Number tested
      num_tested <- tempdata$'NUM_FIRST_VISIT-TOTAL'
      num_tested <- num_tested[num_tested!="A" &  !is.na(num_tested)]
      if(length(num_tested)>=1)
      {
        num_tested <- mean(as.numeric(as.character(num_tested)))
      } else num_tested <- NA
      JaWb_CongenSyphdb_v0$'N tested (1st ANC visit)'[ii] = num_tested

      #Percentage tested
      per_tested <- tempdata$'PER_FIRST_VISIT-TOTAL_Calculated'
      per_tested <- per_tested[per_tested!="A" &  !is.na(per_tested)]
      if(length(per_tested)>=1)
      {
        per_tested <- mean(as.numeric(as.character(per_tested)))
      } else per_tested <- NA
      JaWb_CongenSyphdb_v0$'Syphilis-tested (1st ANC, %)'[ii] = NA #//Screening Sheet

      #Num tested, any visit
      num_tested <- tempdata$'NUM_ANY_VISIT-TOTAL'
      num_tested <- num_tested[num_tested!="A" &  !is.na(num_tested)]
      if(length(num_tested)>=1)
      {
        num_tested <- mean(as.numeric(as.character(num_tested)))
      } else num_tested <- NA
      JaWb_CongenSyphdb_v0$'N tested (any visit ANC)'[ii] = num_tested

      #Num any visits
      num_tested <- tempdata$'DEN_ANY_VISIT-TOTAL'
      num_tested <- num_tested[num_tested!="A" &  !is.na(num_tested)]
      if(length(num_tested)>=1)
      {
        num_tested <- mean(as.numeric(as.character(num_tested)))
      } else num_tested <- NA
      JaWb_CongenSyphdb_v0$'N any ANC visits'[ii] = num_tested

      #Per any visits
      per_tested <- tempdata$'PER_ANY_VISIT-TOTAL_Calculated'
      per_tested <- per_tested[per_tested!="A" &  !is.na(per_tested)]
      if(length(per_tested)>=1)
      {
        per_tested <- mean(as.numeric(as.character(per_tested)))
      } else per_tested <- NA
      JaWb_CongenSyphdb_v0$'Syphilis-tested (any ANC visit,%)'[ii] = per_tested
    }#End if(nrow(tempdata)>=1)
  }


  JaWb_EarlyANC <- data.frame(ISO3 <- JaWbSCreening$COUNTRY_CODE)
  JaWb_EarlyANC$Country <- JaWbSCreening$COUNTRY_NAME
  JaWb_EarlyANC$Region <- JaWbSCreening$REGION_WHO
  JaWb_EarlyANC$'Start year' <- JaWbSCreening$period_start_date
  JaWb_EarlyANC$'End year' <- JaWbSCreening$period_end_date
  JaWb_EarlyANC$'Coverage early ANC' <- JaWbSCreening$Coverage
  JaWb_EarlyANC$'N early ANC' <- JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Code, first ANC threshold time' <- NA#JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Sample size' <- NA#JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Source code' <- NA#JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Source for Source code' <- NA#JaWbSCreening$'NUM_FIRST_VISIT-TOTAL'
  JaWb_EarlyANC$'Early ANC' <- JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Late ANC' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'All women' <- 1#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Before cut-off' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'After cut-off' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'ABO risk, average, treated women' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Time first ANC, before cut-off' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Time first ANC, after cut-off' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Time of first ANC, national average' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100
  JaWb_EarlyANC$'Source' <- NA#JaWbSCreening$'PER_FIRST_VISIT-TOTAL'/100


  JaWb_EarlyANC$SDG_Region <- sapply(JaWb_EarlyANC$ISO3, function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  #########################################################################################################
  #########################################################################################################

  #Reading Eline's data
  ElWb_CongenSyphdb <- openxlsx::read.xlsx(fin_name_cs_db, sheet="CongenSyph db")[,1:30]
  ElWb_CongenSyphdb$COUNTRY[ElWb_CongenSyphdb$COUNTRY=="South-Sudan"] <- "South Sudan"

  names(ElWb_CongenSyphdb) <- c("Rank, 2012 ABO cases", "Rank, 2016 ABO cases", "Rank 2012, CS case RATE", "Rank 2016, CS case RATE", "Country","WHO Region",
                                "ISO code", "ISO3nmb", "Year", "Live Births", "Still births", "Stillbirths, Blencowe & Hogan 2016", "Pregnancies", "Still/Live births",
                                "Source of Live and/or Stillbirths", "Women with >= 1 ANC visit (%)", "Source of ANC1", "N tested (1st ANC visit)", "N, 1st visits",
                                "Syphilis-tested (1st ANC, %)", "N tested (any visit ANC)", "N any ANC visits", "Syphilis-tested (any ANC visit,%)", "Source of Test coverage",
                                "ANC women with syphilis, treated, N", "Syphilis-infected ANC", "Treated (%)", "Source of Treated", "Congenital syphilis case REPORTS",
                                "CS case report rate")

  ElWb_CongenSyphdb$SDG_Region <- sapply(ElWb_CongenSyphdb$'ISO code', function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })


  #Fixing Years
  ElWb_CongenSyphdb <- subset(ElWb_CongenSyphdb,!is.na(Year) & !is.na(ElWb_CongenSyphdb$'ISO code'))
  ElWb_CongenSyphdb$Year[ElWb_CongenSyphdb$Year=="2007-9"] <- "2008"

  ElWb_CongenSyphdb$Year <- as.numeric(as.character(ElWb_CongenSyphdb$Year))
  ElWb_CongenSyphdb$Year[ElWb_CongenSyphdb$Year>3000] <- ElWb_CongenSyphdb$Year[ElWb_CongenSyphdb$Year>=3000]-1000

  #Merging with Jane's data
  ElWb_CongenSyphdb <- rbind(ElWb_CongenSyphdb,JaWb_CongenSyphdb_v0)
  ElWb_CongenSyphdb$imputed_pregnancies <- NA
  ElWb_CongenSyphdb$imputed_ratio_births <- NA
  ElWb_CongenSyphdb$imputed_prevalence <- NA


  Copy_of_ElWb_CongenSyphdb <- ElWb_CongenSyphdb

  ElWb_CongenSyphdb <- ElWb_CongenSyphdb[-c(1:nrow(ElWb_CongenSyphdb)),]

  ######*********************************************************************#####
  #Extend the years
  ######*********************************************************************#####
  ######*
  proj_years <- 1990:2025
  for(tempiso in unique(Copy_of_ElWb_CongenSyphdb$`ISO code`))
  {
    temp_ctr <- subset(Copy_of_ElWb_CongenSyphdb,Copy_of_ElWb_CongenSyphdb$`ISO code`==tempiso)
    if(nrow(temp_ctr)==0) next
    missing_years <- proj_years[which(!is.element(proj_years,temp_ctr$Year))]
    if(length(missing_years)>=1)
    {
      temp_ctr_add <- do.call("rbind", replicate(length(missing_years), temp_ctr[1,], simplify = FALSE))
      temp_ctr_add[,] <- NA
      temp_ctr_add$Year <- missing_years
      temp_ctr_add$Country <- temp_ctr$Country[1]
      temp_ctr_add$`WHO Region` = temp_ctr$`WHO Region`[1]
      temp_ctr_add$`ISO code` <- temp_ctr$`ISO code`[1]
      temp_ctr_add$ISO3nmb <- temp_ctr$ISO3nmb[1]
      temp_ctr <- rbind(temp_ctr,temp_ctr_add)
      rm(temp_ctr_add)
    }

    temp_ctr <- subset(temp_ctr,!is.na(Year))
    idx_sorted_years <- sort.int(temp_ctr$Year,index.return = T)$ix
    temp_ctr <- temp_ctr[idx_sorted_years,]

    temp_ctr$imputed_pregnancies <- ifelse(is.na(temp_ctr$Pregnancies),TRUE,FALSE)
    temp_ctr$imputed_ratio_births <- ifelse(is.na(temp_ctr$imputed_ratio_births),TRUE,FALSE)

    ElWb_CongenSyphdb <- rbind(ElWb_CongenSyphdb,temp_ctr)
  }



  ElWb_CongenSyphdb$'Still births' <- sapply(ElWb_CongenSyphdb$'Still births', function(x)
  {
    res <- NA;
    if(!is.na(x)) res <- as.numeric(as.character(x))
    res
  })

  ElWb_CongenSyphdb$'Live Births' <- sapply(ElWb_CongenSyphdb$'Live Births', function(x)
  {
    res <- NA;
    if(!is.na(x)) res <- as.numeric(as.character(x))
    res
  })

  ElWb_CongenSyphdb$'Still/Live births'[] <- NA;
  ElWb_CongenSyphdb$'Still/Live births' <- ElWb_CongenSyphdb$'Still births'/ElWb_CongenSyphdb$'Live Births'

  #options(warn=2)
  #
  ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)' <- sapply(ElWb_CongenSyphdb$'Women with >= 1 ANC visit (%)', function(x)
  {
    res <- NA;
    if(!is.na(x))
    {
      cpx <- x
      while(grepl("or",cpx,TRUE))
      {
        idx <- unlist(gregexpr("or",cpx))[1]
        cpx <- substr(cpx,idx+2,nchar(cpx))
      }

      while(grepl("\\?",cpx,TRUE))
      {
        idx <- unlist(gregexpr("\\?",cpx))[1]
        if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
      }
      #cat(x, ", cpx=",cpx,"\n")
      while(grepl("-",cpx,TRUE))
      {
        idx <- unlist(gregexpr("-",cpx))[1]
        if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
      }
      while(grepl("%",cpx,TRUE))
      {
        idx <- unlist(gregexpr("%",cpx))[1]
        if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
      }
      res <- as.numeric(as.character(cpx))
    }
    res
  })

  gm_convertchartonum <- function(cpx)
  {
    while(grepl("or",cpx,TRUE))
    {
      idx <- unlist(gregexpr("or",cpx))[1]
      cpx <- substr(cpx,idx+2,nchar(cpx))
    }

    while(grepl("\\?",cpx,TRUE))
    {
      idx <- unlist(gregexpr("\\?",cpx))[1]
      if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
    }
    #cat(x, ", cpx=",cpx,"\n")
    while(grepl("-",cpx,TRUE))
    {
      idx <- unlist(gregexpr("-",cpx))[1]
      if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
    }
    while(grepl("%",cpx,TRUE))
    {
      idx <- unlist(gregexpr("%",cpx))[1]
      if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
    }

    for(tch in c(">","<",">=", "<=", "="))

      while(grepl(">",cpx,TRUE))
      {
        idx <- unlist(gregexpr(">",cpx))[1]
        if(idx>1) cpx <- substr(cpx,1,nchar(cpx)-1) else if(idx==1) cpx <- substr(cpx,2,nchar(cpx))
      }
    res <- as.numeric(as.character(cpx))
  }

  ####################################
  ###Syphilis-tested (1st ANC, %)
  ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)' <- sapply(ElWb_CongenSyphdb$'Syphilis-tested (1st ANC, %)', function(x)
  {
    res <- NA;
    if(!is.na(x))
    {
      res <- gm_convertchartonum(x)#as.numeric(as.character(x))
    }
    res
  })


  #####################################
  ###`Treated (%)`
  ElWb_CongenSyphdb$'Treated (%)' <- sapply(ElWb_CongenSyphdb$'Treated (%)', function(x)
  {
    res <- NA;
    if(!is.na(x)& (x!="No positive cases to be treated."))
    {
      res <- gm_convertchartonum(x)#as.numeric(as.character(x))
    }
    res
  })

  #Variable to be converted to numeric
  listvar <- c("Rank, 2012 ABO cases", "Rank, 2016 ABO cases", "Rank 2012, CS case RATE", "Rank 2016, CS case RATE",  #"Country", "WHO Region", "ISO code", "ISO3nmb",
               "Year", "Live Births", "Still births", "Stillbirths, Blencowe & Hogan 2016",
               "Pregnancies", "Still/Live births",  "Women with >= 1 ANC visit (%)", #"Source of ANC1", "Source of Live and/or Stillbirths",
               "N tested (1st ANC visit)", "N, 1st visits", "Syphilis-tested (1st ANC, %)",  "N tested (any visit ANC)", "N any ANC visits",
               "Syphilis-tested (any ANC visit,%)",  "ANC women with syphilis, treated, N", "Syphilis-infected ANC", #"Source of Test coverage",
               "Treated (%)",  "Congenital syphilis case REPORTS", "CS case report rate")#, "Source of Treated",

  for(vname in listvar)
  {
    jj <- which(names(ElWb_CongenSyphdb)==vname)
    ElWb_CongenSyphdb[,jj] <- sapply(ElWb_CongenSyphdb[,jj], function(x)
    {
      res <- NA;
      if(!is.na(x)& (x!="No positive cases to be treated.") & (!grepl("11 of 12 clinics",x)) & (!grepl("Near-universal, near-100",x))&
         (!grepl("https://www.rivm.nl",x))&(!grepl("Numerator is sentinel survey",x))&(!grepl("83%, 64%, 77%, 69%, 76%" ,x))&
         (!grepl("Wijesooriya-S et al." ,x))&(!grepl("Wilkinson-D TMIH" ,x)))
      {
        if(x=="FALSE" | x=="false")
        {
          res <- "0"
        } else if(x=="TRUE"| x=="true")
        {
          res <- "0"
        } else
        {
          res <- gm_convertchartonum(x)#as.numeric(as.character(x))
        }
      }
      res
    })
  }

  ################################################################################
  ######Imputation of pregnancies#################################################
  ################################################################################


  SpectrumBirths <- openxlsx::read.xlsx(fin_name_allbirths, startRow = 2)
  #Countries with missing Spectrum files, we interpolate using a linear model
  for(ctr in unique(ElWb_CongenSyphdb$Country))
  {
    idx_ctr <- which(ElWb_CongenSyphdb$Country==ctr)
    if(any(is.na(ElWb_CongenSyphdb$Pregnancies[idx_ctr])))
    {
      idx_na <-  which(is.na(ElWb_CongenSyphdb$Pregnancies) & ElWb_CongenSyphdb$Country==ctr & !is.na(ElWb_CongenSyphdb$Year))
      idx_notna <-  which(!is.na(ElWb_CongenSyphdb$Pregnancies) & ElWb_CongenSyphdb$Country==ctr & !is.na(ElWb_CongenSyphdb$Year))
      if(length(idx_notna)>=2 & length(idx_na)>=1)
      {
        mod <- lm(Pregnancies~Year,data=ElWb_CongenSyphdb[idx_notna,])
        pred <- predict(mod, newdata=ElWb_CongenSyphdb[idx_na,])
        pred[pred<0] <- 0
        ElWb_CongenSyphdb$Pregnancies[idx_na] <- pred
      }
    }#End if(any(is.na(ElWb_CongenSyphdb$Pregnancies[idx_ctr])))
  }


  for(ii in which(is.na(ElWb_CongenSyphdb$Pregnancies) | is.na(ElWb_CongenSyphdb$`Live Births`)))
  {
    iso = ElWb_CongenSyphdb$`ISO code`[ii]
    year_ii <- ElWb_CongenSyphdb$Year[ii]
    if(is.na(year_ii)) next
    ratiostilllive <- ElWb_CongenSyphdb$`Still/Live births`[ii]
    livebirths <- SpectrumBirths[SpectrumBirths$ISO3==iso, names(SpectrumBirths) == as.character(year_ii)]
    if(length(livebirths)>=1)
    {
      livebirths <- sum(livebirths)
      still <- ratiostilllive/(1+ratiostilllive)*livebirths;
      live <- livebirths/(1+ratiostilllive);
      ElWb_CongenSyphdb$`Live Births`[ii] <- live;
      ElWb_CongenSyphdb$`Still births`[ii] <- still;
      ElWb_CongenSyphdb$Pregnancies[ii] <- still+live;
    } else
    {
      livebirths <- ElWb_CongenSyphdb$Pregnancies[ii]
      still <- ratiostilllive/(1+ratiostilllive)*livebirths;
      live <- livebirths/(1+ratiostilllive);
      ElWb_CongenSyphdb$`Live Births`[ii] <- live;
      ElWb_CongenSyphdb$`Still births`[ii] <- still;
    }
  }

  ElWb_CongenSyphdb$SDG_Region <- sapply(ElWb_CongenSyphdb$'ISO code', function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })
  ################################################################################
  ################################################################################
  ################################################################################


  ElWb_RiskABO_Asumptions <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Assumptions",rows=2:5, cols=1:10)
  names(ElWb_RiskABO_Asumptions)<- c("X1","treated mothers", "CS, risk: liveborn", "CS risk: stillbirth", "CS risk: neonatal death",
                                     "CS risk: prematurity or LBW", "All outcomes", "LB on ABO risk assumption", "UB on ABO risk assumption",
                                     "Variance on ABO risk assumption")

  ElWb_SE_Asumptions <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Assumptions",rows=11:17, cols=1:3)
  ElWb_ABO_Asumptions <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Assumptions",rows=10:17, cols=10:14)
  names(ElWb_ABO_Asumptions) <- c("WHO Regions", "Indicator","ABO risk, treated mothers","Average timing of first ANC", "Number of Pregnancies, 2016")


  #
  ElWb_TimingANC_Asumptions <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Assumptions",rows=20:30, cols=1:10)
  names(ElWb_TimingANC_Asumptions) <- c("Code, in data", "Timing of first ANC", "Median weeks in the pregnancy, of first ANC", "ABO risk probability/treated syphilis-infected mother",
                                        "Minimum (weeks)", "Maximum (weeks)", "Average (weeks)", "ABO risk probability/untreated syphilis-infected mother", "Added risk, from early to late",
                                        "Frequency")

  ElWb_EarlyANC <- openxlsx::read.xlsx(fin_name_cs_db,sheet="Early ANC",rows=1:528, cols=1:21)
  names(ElWb_EarlyANC) <- c("ISO3", "Country", "Region", "Start year", "End year", "Coverage early ANC", "N early ANC", "Code, first ANC threshold time",
                            "Sample size", "Source code", "Source for Source code", "Early ANC", "Late ANC", "All women", "Before cut-off", "After cut-off",
                            "ABO risk, average, treated women", "Time first ANC, before cut-off", "Time first ANC, after cut-off",
                            "Time of first ANC, national average", "Source")


  ElWb_EarlyANC$SDG_Region <- sapply(ElWb_EarlyANC$'ISO3', function(x)
  {
    res <- NA
    idx <- which(SDGRegions$COUNTRY_CODE==x)
    if(length(idx)>=1) res <- SDGRegions$'SDG Regions'[idx[1]]
    res
  })

  ElWb_EarlyANC$`Start year`[is.na(ElWb_EarlyANC$`Start year`)] <- ElWb_EarlyANC$`End year`[is.na(ElWb_EarlyANC$`Start year`)]
  ElWb_EarlyANC <- subset(ElWb_EarlyANC, (!is.na(ElWb_EarlyANC$`Start year`)) & (!is.na(ElWb_EarlyANC$`End year`)))


  LongEarlyANCData <- data.frame();
  for(ctriso in unique(ElWb_EarlyANC$ISO3))
  {
    if(is.na(ctriso)) next;
    long_ctr_data <- data.frame()
    ctr_data <- subset(ElWb_EarlyANC,ISO3==ctriso)
    for(ii in 1:nrow(ctr_data))
    {
      ctr_dat_slice <- ctr_data[ii,]
      n_years <- ctr_dat_slice$`End year`[1]-ctr_dat_slice$`Start year`[1]+1
      if(n_years<0)
      {
        n_years = -n_years;
        fy <- ctr_dat_slice$`Start year`[1]
        ctr_dat_slice$`Start year`[1] =ctr_dat_slice$`End year`[1]
        ctr_dat_slice$`End year`[1] = fy
      }
      temp_long_ctr_data <- do.call("rbind", replicate(n_years, ctr_dat_slice, simplify = FALSE))
      temp_long_ctr_data$`Start year` = ctr_dat_slice$`Start year`[1]+0:(n_years-1)
      temp_long_ctr_data$`End year` = ctr_dat_slice$`Start year`[1]+1:n_years;
      long_ctr_data <- rbind(long_ctr_data,temp_long_ctr_data)
    }

    nanames <- c("Coverage early ANC","N early ANC","Code, first ANC threshold time", "Sample size",
                 "Source code","Source for Source code","Early ANC","Late ANC","All women","Before cut-off",
                 "After cut-off","ABO risk, average, treated women","Time first ANC, before cut-off","Time first ANC, after cut-off",
                 "Time of first ANC, national average","Source","SDG_Region")

    temp_proj_years <- 1980:2025
    fin_long_ctr_data <- data.frame();
    for(yy in temp_proj_years)
    {
      temp_long_ctr_data <- subset(long_ctr_data,long_ctr_data$'Start year'==yy)
      if(nrow(temp_long_ctr_data)==0)
      {
        temp_long_ctr_data <- long_ctr_data[1,]
        temp_long_ctr_data$`Start year` = yy
        temp_long_ctr_data$`End year` = yy+1
        temp_long_ctr_data[,is.element(names(temp_long_ctr_data),nanames)] <- NA
      } else if (nrow(temp_long_ctr_data)>1)
      {
        temp_long_ctr_data0 <- temp_long_ctr_data
        temp_long_ctr_data <- temp_long_ctr_data[1,]
        for(vname in nanames[c(1,2,4,7,8,9,10,11,12,13,14,15)])
        {
          temp_long_ctr_data[,names(temp_long_ctr_data)==vname] <- mean(temp_long_ctr_data0[,names(temp_long_ctr_data0)==vname],na.rm=T)
        }
      }
      fin_long_ctr_data <- rbind(fin_long_ctr_data,temp_long_ctr_data)
    }

    LongEarlyANCData <- rbind(LongEarlyANCData,fin_long_ctr_data)
  }


  CongenDataRaw <- ElWb_CongenSyphdb; rm(ElWb_CongenSyphdb)
  CongenDataRaw <- subset(CongenDataRaw, !(Country%in%c("AFR total","AMR total","EMR total","EUR total","SEAR total","WPR total")))
  CongenDataRaw$Country[CongenDataRaw$Country=="Saint Vincent and Grenadines"] = "Saint Vincent and the Grenadines"
  CongenDataRaw$Country[CongenDataRaw$Country=="Bahamas (the)"] = "Bahamas"
  CongenDataRaw$Country[CongenDataRaw$Country=="Cabo Verde"] = "Cape Verde"
  CongenDataRaw$Country[CongenDataRaw$Country=="Central African Republic (the)"] = "Central African Republic"
  CongenDataRaw$Country[CongenDataRaw$Country=="Comoros (the)"] = "Comoros"
  cst <- CongenDataRaw$Country[grep("ville",CongenDataRaw$Country)]
  if(length(cst)>=1)
  {
    CongenDataRaw$Country[CongenDataRaw$Country==cst[1]] = "Congo"
  }
  CongenDataRaw$Country[CongenDataRaw$Country=="Congo (Brazzaville)"] = "Congo"
  CongenDataRaw$Country[CongenDataRaw$Country=="C?te d'Ivoire"] = "Cote d'Ivoire"
  CongenDataRaw$Country[CongenDataRaw$Country=="Dem. People's Republic of Korea"] = "Democratic People's Republic of Korea"
  CongenDataRaw$Country[CongenDataRaw$Country=="Dominican Republic (the)"] = "Dominican Republic"
  CongenDataRaw$Country[CongenDataRaw$Country=="Swaziland"] = "Eswatini"
  CongenDataRaw$Country[CongenDataRaw$Country=="Micronesia (Fed. States of)"] = "Micronesia (Federated States of)"
  CongenDataRaw$Country[CongenDataRaw$Country=="The Gambia"] = "Gambia"
  CongenDataRaw$Country[CongenDataRaw$Country=="United Kingdom of Great Britain and Northern Ireland"] = "United Kingdom"
  CongenDataRaw$Country[CongenDataRaw$Country=="Venezuela, Bolivarian Republic of"] = "Venezuela (Bolivarian Republic of)"
  CongenDataRaw$Country[CongenDataRaw$Country=="Viet Nam"] = "Vietnam"
  CongenDataRaw$Country[CongenDataRaw$Country=="Czech Republic"] = "Czechia"


  result <- list(CongenDataRaw=CongenDataRaw,SyphDataRaw=SyphDataRaw)
  return(result)
}
