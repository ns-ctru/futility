#' Conditional Power Calculation for Futility
#'
#' @description Re-calculate power based on observed summary statistics.
#'
#' @details
#'
#' @param df Data frame of recruited data, must contain the \code{outcome}, \code{group} and optionally \code{recruited}.
#' @param outcome The otucome variable you are assessing.
#' @param group Binary indicator of allocation.
#' @param recruited Date of an individuals recruitment.
#' @param planned.case Number of planned cases to be recruited.
#' @param planned.control Number of planned controls to be recruited.
#' @param n.recruited.case Number of recruited cases to date.
#' @param n.recruited.control Number of recruited controls to date.
#' @param mean.recruited.case Mean of outcome in cases recruited to date.
#' @param mean.recruited.control Mean of outcome in controls recruited to date.
#' @param sd.recruited.case Standard deviation in cases recruited to date.
#' @param sd.recruited.control Standard deviation in controls recruited to date.
#' @param sd.recruited.pooled Standard deviation in cases and controls recruited to date.
#' @param method Type of outcome \code{continuous} (default) | \code{binary} | \code{survival}
#' @param assumed.delta Assumed treatment effect.
#' @param assumed.delta.final Assumed treatment effect at the end of the trial (i.e. the estimate used in the original sample size calculation).
#' @param futility.threshold Pre-specified conditional power/futility threshold for deciding to stop.
#' @param alpha Planned Type I error rate.
#' @param beta Planned Type II error rate.
#' @param null.treatment Difference between groups if the null is true.
#'
#' @return A list of results depending on the options specified.
#'
#' @examples
#'
#'
#'
#' @references
#'
#'  Lachin JM (2005) A review of methods for futility stopping based on conditional power. Statistics in Medicine 24(18):2747-2764
#'  Proschan MA, Lan KKG, and Wittes JT (2006) Statistical monitoring of clinical trials: a unified approach. Springer Science & Business Media
#' @export
futility <- function(df                     = data,
                     outcome                = value,
                     group                  = rand,
                     recruited              = recruit.date,
                     n.planned              = c(100, 100),
                     method                 = 'continuous',
                     assumed.delta          = 1,
                     assumed.delta.final    = 1,
                     futility.threshold     = NA,
                     alpha                  = 0.05,
                     beta                   = 0.9,
                     null.treatment         = 0,
                     ...){
    ## ## Check options
    ## if(is.na(planned.case) | is.null(planned.case) | is.na(planned.control) | is.null(planned.control)){
    ##     print('Error : You must specify the planned number of participants for both arms.')
    ##     exit
    ## }
    ## if(is.na(n.recruited.case) | is.null(n.recruited.case) | is.na(n.recruited.control) | is.null(n.recruited.control)){
    ##     print('Error : You must specify the number of participants already recruited for both arms.')
    ##     exit
    ## }
    ## if(is.na(recruited.case) | is.null(recruited.case) | is.na(recruited.control) | is.null(recruited.control)){
    ##     print('Error : This is the immediate version of the command and you must specify the mean/proportion in both arms.  If you wish to use an existing data set and have these calculated automatically use the futility_data() function.')
    ##     exit
    ## }
    ## if(!is.integer(planned.case)){
    ##     print('Error : The planned number of cases should be an integer.')
    ## }
    ## if(!is.integer(planned.control)){
    ##     print('Error : The planned number of controls should be an integer.')
    ## }
    ## if(!is.numeric(alpha) | !is.numeric(beta)){
    ##     print('Error : alpha and beta must be numeric.')
    ## }
    ## if(alpha <0 | alpha > 1){
    ##     print('Error : alpha must be a numerical value (0 <= alpha <= 1)')
    ## }
    ## if(beta <0 | beta > 1){
    ##     print('Error : beta must be a numerical value (0 <= beta <= 1)')
    ## }
    ## ## Check recruited values
    ## if(method == 'binary'){
    ##     if((recruited.case < 0 | recruited.case > 1)){
    ##         print('Error : The proportion in the case arm must be 0 <= recruited.case <= 1 ')
    ##         exit
    ##     }
    ##     if((recruited.control < 0 | recruited.control > 1)){
    ##         print('Error : The proportion in the control arm must be 0 <= recruited.control <= 1 ')
    ##         exit
    ##     }
    ## }
    ## Initialise list of results for returning
    results <- list()
    #############################################################################
    ## CONTINUOUS OUTCOMES                                                     ##
    #############################################################################
    if(method == 'continuous'){
        ## Summarise data at current interim
        if(!is.null(df) & !is.null(group) & !is.null(outcome)){
            summary.time <- group_by(df, group) %>%
                            summarise(mean = mean(outcome, na,rm = TRUE),
                                      sd   = sd(outcome, na.rm = TRUE),
                                      n    = n())
            ## Derive SE
            summary.time <- mutate(summary.time,
                                   se = sd / sqrt(n))
        }
        ## Cumulative summary of data (only possible if recruited date is provided)
        if(!is.null(df) & !is.null(group) & !is.null(outcome) & !is.null(recruited)){
            ## TODO - How to group by sequential dates and summarise events upto that time point?
            ##        cummean() provides some functionality, but no SD
            ##        https://stackoverflow.com/questions/34874582/calculating-cumulative-standard-deviation-by-group-using-r
            summary.over.time <- group_by(df, group) %>%
                                 summarise(mean = cummean(outcome, na,rm = TRUE),
                                           sd   = sd(outcome, na.rm = TRUE),
                                           n    = n())
        }
        ## Expected Theta = z-score - drift under H0
        results$theta         <- qnorm(1 - (alpha / 2)) + qnrom(1 - beta)
        results$interim.h0    <- null.treatment
        ## Interim treatment effect and standard deviation
        results$interim.mean.diff <- summary.time$mean[1] - summary.time$mean[2]
        results$interim.sd.diff   <- sqrt((summary.time$sd[1]^2 / summary.time$n[1]) + (summary.time$sd[2]^2 / summary.time$n[2]))
        ## Standardised treatment effect
        results$interim.z <- results$interim.mean.diff / results$interim.sd.diff
        #############################################################################
        ## Compute Genearlised information Fraction (Lachin (2005) eq1)            ##
        #############################################################################
        results$information.fraction <- ((summary.time$n[1]^-1 + summary.time$n[2]^-1)^-1) / ((planned.case^-1 + planned.control^-1)^-1)
        ## if(results$information.fraction < 0 | results$information.fraction > 1){
        ##     print('Error : Something has gone wrong, the information fraction should be in the range of 0 to 1.  Please check you have correctly specified the planned and recruited numbers.')
        ##     exit
        ## }
    }
    #############################################################################
    ## BINARY OUTCOMES                                                         ##
    #############################################################################
    ## TODO - Once Continuous outcomes sorted
    else if(method == 'binary'){
    }
    #############################################################################
    ## SURVIVAL OUTCOMES                                                       ##
    #############################################################################
    ## TODO - Once Continuous outcomes sorted
    else if(method == 'survival'){
    }
    ## Return results
    return(results)
}
