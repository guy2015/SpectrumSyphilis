#' Outputs default Diagnostic Tests parameters
#' @return A data frame with Adjusment factors for diagnostic tests.
#' @examples get_DefaultDiagTests()
get_DefaultDiagTests <- function()
{
   return(Default_DiagnosticTest)
}

#' Outputs default Syphilis durations
#' @return A list containing syphilis durations by WHO regions.
#' @examples get_DefaultDiagTests()
get_DefaultSyphDur <- function()
{
  return(list(Sypdur_A = 1.28, Sypdur_B = 2.42, Sypdur_C = 4.13))
}

#' Outputs default Prevalence Ratios to be used when prevalence among KPs are unavailable
#' @return A list containing syphilis durations by WHO regions.
#' @examples get_DefaultPOR()
get_DefaultPOR <- function()
{
  result <- list(MtoMSM = list(logPOR=0.09984533, sdlogPOR=0.1), FtoFSW=list(logPOR=0.09984533, sdlogPOR=0.1))
  return(result)
}
