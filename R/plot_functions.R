#' Dot plots for influence diagnostics
#'
#' This is a function that can be used to create (modified) dotplots for the
#' diagnostic measures.  The plot allows the user to understand the distribution
#' of the diagnostic measure and visually identify unusual cases.
#' 
#' @note
#' The resulting plot uses \code{coord_flip} to rotate the plot, so when
#' adding customized axis labels you will need to flip the names 
#' between the x and y axes. 
#'
#' @param x values of the diagnostic of interest
#' @param cutoff value(s) specifying the boundary for unusual values of the 
#' diagnostic. The cutoff(s) can either be supplied by the user, or automatically
#' calculated using measures of internal scaling if \code{cutoff = "internal"}.
#' @param name what diagnostic is being plotted
#' (one of \code{"cooks.distance"}, \code{"mdffits"}, \code{"covratio"}, 
#' \code{"covtrace"}, \code{"rvc"}, or \code{"leverage"}).
#' This is used for the calculation of "internal" cutoffs.
#' @param data data frame to use (optional)
#' @param index optional parameter to specify index (IDs) of \code{x} values.
#' If \code{NULL}(default), values will be indexed in the order of the vector 
#' passed to \code{x}. 
#' @param modify specifies the \code{geom} to be used to produce a 
#' space-saving modification: either \code{"dotplot"} or \code{"boxplot"}
#' @param ... other arguments to be passed to \code{ggplot()}
#' @author Adam Loy \email{loyad01@@gmail.com}
#' @examples 
#' data(sleepstudy, package = 'lme4')
#' fm <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' 
#' #Observation level deletion and diagnostics
#' obs.infl <- hlm_influence(fm, level = 1)
#' 
#' dotplot_diag(x = obs.infl$cooksd, cutoff = "internal", name = "cooks.distance", modify = FALSE)
#' 
#' dotplot_diag(x = obs.infl$mdffits, cutoff = "internal", name = "cooks.distance", modify = FALSE)
#'
#' # Subject level deletion and diagnostics
#' subject.infl  <- hlm_influence(fm, level = "Subject")
#' 
#' dotplot_diag(x = subject.infl$cooksd, cutoff = "internal",
#'              name = "cooks.distance", modify = FALSE)
#'              
#'dotplot_diag(x = subject.infl$mdffits, cutoff = "internal", name = "mdffits", modify = "dotplot")
#' @keywords hplot
#' @importFrom rlang .data
#' @export
dotplot_diag <- function(x, cutoff, 
                         name = c("cooks.distance", "mdffits", "covratio", "covtrace", "rvc", "leverage"),
                         data, index = NULL, modify = FALSE, ...) {
  
  # Address no visible binding for R CMD CHECK 
  n.factor <- n <- value <- extreme <- NULL
  
  if(!modify %in% c(FALSE, "boxplot", "dotplot")) {
    stop("modify should be FALSE or either 'boxplot' or 'dotplot'")
  }
  
  if(modify != FALSE & missing(cutoff)){
    stop("a cutoff should be specified if a modified dotplot is requested")
  }
  
  if(modify != FALSE & missing(name)){
    stop("a name should be specified if a modified dotplot is requested")
  }
  
  if(!missing(cutoff)){
    if( !is.numeric(cutoff) && cutoff != "internal" ){
      stop("cutoff should be numeric or 'internal'")
    }
  } else{
    cutoff <- NULL
  }
  
  if(!missing(name)){
    if(!name %in% c("cooks.distance", "mdffits", "covratio",
                    "covtrace", "rvc", "leverage")) {
      stop("name should be one of 'cooks.distance', 'mdffits', 'covratio', 
           'covtrace', 'rvc', 'leverage'")
    }
    }
  
  if(missing(data)) {
    data <- data.frame() 
    if(is.null(index)) index <- factor(seq(1, length(x)))
  } else{
    if(is.null(index)) index <- factor(seq(1, nrow(data)))
  }
  
  name <- match.arg(name)
  x <- eval(substitute(x), data, parent.frame())
  index <- eval(substitute(index), data, parent.frame())
  
  if(class(x) %in% c("fixef.dd", "vcov.dd")) x <- as.numeric(x)

  df <- data.frame(index = index, value = x) 
  df <- df %>% 
    arrange(desc(.data$value)) %>%
    mutate(n = nrow(df):1, n.factor = factor(seq(nrow(df), 1))) 
  
  if(!is.null(cutoff)){
    
    if(!is.numeric(cutoff)) cutoff <- internal_cutoff(x = x, name = name)
    
    if(is.numeric(cutoff)){ 
      if(!name %in% c("covratio", "rvc") ){
        df <- df %>%
          mutate(extreme = ifelse(.data$value > cutoff, TRUE, FALSE))
      }
      else {
        df <- df %>%
          mutate(extreme = ifelse(.data$value < cutoff[1] | .data$value > cutoff[2], TRUE, FALSE))
      }
      
      
      if (modify == FALSE) {
        p <- ggplot(dplyr::filter(df, .data$extreme == FALSE)) + 
          geom_point(aes(x = n, y = value), colour = I("blue"), inherit.aes = FALSE) 
        
        if( sum(df$extreme) > 0 ){
          p <- p + 
            geom_point(
              data = subset(df, extreme == TRUE), 
              aes(x = n, y = value), 
              color = I("red"), 
              shape = 17
            ) +
            ggrepel::geom_text_repel(
              data = df[1:5,], 
              aes(x = n, y = value, label = index,  hjust=.5, vjust=1.5, size=3)
            ) 
        }
        
        if (!name %in% c("covratio", "rvc")) {
          p + geom_hline(aes(yintercept = cutoff), colour = I("red")) +
            theme(
              legend.position = "none", 
              axis.title.y = element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank()
            ) +
            labs(y = name) +
            coord_flip() 
        }
        else {
          p + geom_hline(aes(yintercept = cutoff[1]), colour = I("red")) +
            geom_hline(aes(yintercept = cutoff[2]), colour = I("red")) +
            theme(
              legend.position = "none", 
              axis.title.y = element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank()
            ) +
            labs(y = name) +
            coord_flip() 
        }
        
      }
      
      else { #space saving measure 
        levels(df$n.factor)[levels(df$n.factor) %in% filter(df, .data$extreme == FALSE)$n.factor] <- 0
        
        p <- ggplot(dplyr::filter(df, .data$extreme == FALSE)) 
        
        if (modify == "boxplot") {
          p <- p + geom_boxplot(aes(x = n.factor, y = value, group = 1), inherit.aes = FALSE) 
        }
        if (modify == "dotplot") {
          p <- p  + geom_point(
            aes(x = n.factor, y = value), 
            colour = I("blue"), 
            inherit.aes = FALSE
          ) 
        }
        
        if( sum(df$extreme) > 0 ){
          
          p <- p +
            geom_point(
              data = filter(df, extreme == TRUE), 
              aes(x = n.factor, y = value),
              color = I("red"), 
              shape = 17, 
              inherit.aes = FALSE
            ) +
            ggrepel::geom_text_repel(
              data = filter(df, .data$extreme == TRUE)[1:5,], 
              aes(x = n.factor, y = value, label = index,  
                  hjust = .5, vjust = 1.5, size = 3), 
              inherit.aes = FALSE
            )  
        }
        
        
        if (!name %in% c("covratio", "rvc")) {
          p + geom_hline(aes(yintercept = cutoff), colour=I("red")) +
            theme(
              legend.position = "none", 
              axis.title.y = element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank()
            ) +
            labs(y = name) +
            coord_flip() 
        }
        else {
          p + geom_hline(aes(yintercept = cutoff[1]), colour=I("red")) +
            geom_hline(aes(yintercept = cutoff[2]), colour=I("red")) +
            theme(
              legend.position = "none", 
              axis.title.y = element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank()
            ) +
            labs(y = name) +
            coord_flip() 
        }
      }  
    }
  }
  else {
    ggplot(df) + geom_point(aes(x = n, y = value)) + coord_flip()
  }
  }



# Calculating a cutoff value for diagnostic measures
#
# This function provides cutoff values using internal scaling. 
# In other words, a measure of relative standing is used (3*IQR)
# to specify unusual values of the diagnostic of interest relative
# to the vector of diagnostics.
#
# @param x a vector
# @param name specification of which diagnostic to plot
internal_cutoff <- function(x, name){
  q3 <- quantile(x, p=0.75)
  x.iqr <- IQR(x)
	
  if(name %in% c("covratio", "rvc")){
    q1 <- quantile(x, p=0.25)
    cutoff <- c(lower = q1 - 3 * x.iqr, upper = q3 + 3 * x.iqr)	
  }
	
  else{cutoff <- q3 + 3 * x.iqr}
	
  return(cutoff)	
}
