#' Conditional Power Calculation for Futility
#'
#' @description Re-calculate power based on observed summary statistics.
#'
#' @details
#'
#' @param planned.case
#' @param planned.control
#' @param recruited.case
#' @param recruited.control
#' @param planned.case
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
futility <- function(planned.case         = 100,
                     planned.control      = 100,
                     n.recruited.case     = 50,
                     n.recruited.control  = 50,
                     recruited.case       = 1,
                     recruited.control    = 2,
                     sd.recruited.case    = 1,
                     sd.recruited.control = 1,
                     sd.recruited.pooled  = 1,
                     method               = 'binary',
                     assumed.delta        = 1,
                     assumed.delta.final  = 1,
                     sd.pooled            = 1,
                     futility.theshold    = NA,
                     alpha                = 0.05,
                     beta                 = 0.9,
                     ...){
    ## Check options
    if(is.na(planned.case) | is.null(planned.case) | is.na(planned.control) | is.null(planned.control)){
        print('Error : You must specify the planned number of participants for both arms.')
        exit
    }
    if(is.na(n.recruited.case) | is.null(n.recruited.case) | is.na(n.recruited.control) | is.null(n.recruited.control)){
        print('Error : You must specify the number of participants already recruited for both arms.')
        exit
    }
    if(is.na(recruited.case) | is.null(recruited.case) | is.na(recruited.control) | is.null(recruited.control)){
        print('Error : This is the immediate version of the command and you must specify the mean/proportion in both arms.  If you wish to use an existing data set and have these calculated automatically use the futility_data() function.')
        exit
    }
    if(!is.integer(planned.case)){
        print('Error : The planned number of cases should be an integer.')
    }
    if(!is.integer(planned.control)){
        print('Error : The planned number of controls should be an integer.')
    }
    if(!is.numeric(alpha) | !is.numeric(beta)){
        print('Error : alpha and beta must be numeric.')
    }
    if(alpha <0 | alpha > 1){
        print('Error : alpha must be a numerical value (0 <= alpha <= 1)')
    }
    if(beta <0 | beta > 1){
        print('Error : beta must be a numerical value (0 <= beta <= 1)')
    }



    }
    ## Check recruited values
    if(method == 'binary'){
        if((recruited.case < 0 | recruited.case > 1)){
            print('Error : The proportion in the case arm must be 0 <= recruited.case <= 1 ')
            exit
        }
        if((recruited.control < 0 | recruited.control > 1)){
            print('Error : The proportion in the control arm must be 0 <= recruited.control <= 1 ')
            exit
        }
    }
    ## Initialise list of results for returning
    results <- list()
    #############################################################################
    ## Compute Genearlised information Fraction (Lachin (2005) eq1)            ##
    #############################################################################
    results$information.fraction <- ((n.recruited.case^-1 + n.recruited.control^-1)^1) / ((planned.case^-1 + planned.control^-1)^-1)
    if(results$information.fraction < 0 | results$information.fraction > 1){
        print('Error : Something has gone wrong, the information fraction should be in the range of 0 to 1.  Please check you have correctly specified the planned and recruited numbers.')
        exit
    }
    #############################################################################
    ## Expected Theta = z-score - drift under H0                               ##
    #############################################################################
    results$theta <- inverse.gaussian(1 - (alpha / 2)) + inverse.gaussian(1 = beta)
    #############################################################################
    ## Interim Treatment effect                                                ##
    #############################################################################
    ## ToDo - query logic of dividing by zero?
    ## results$interim.effect <- (recruited.case - recruited.control) /
    ## Return results
    return(results)
}
