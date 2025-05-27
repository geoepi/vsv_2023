phylodynamic_process <- function(tree, stats, x_limits = c("2010-01-01", "2017-06-01")){
  
  root_age <- stats %>%
    filter(Parameter == "age.root.") %>%
    select(Median) %>%
    pull()
  
  tree_mrsd <- stats %>%
    filter(Parameter == "treeModel.rootHeight") %>%
    select(Median) %>%
    pull() + root_age
  
  # tree$root.time <- root_age 
  
  coal_pref <- phylodyn::BNPR_PS(tree, lengthout = 500, 
                                 prec_alpha = 0.001, 
                                 prec_beta = 0.001,
                                 beta1_prec = 0.0001, 
                                 fns = NULL, 
                                 log_fns = FALSE, 
                                 simplify = TRUE,
                                 derivative = FALSE, 
                                 forward = TRUE)
  
  coal_pref_df <- as.data.frame(
    cbind(
      date = coal_pref$x,
      Ne = coal_pref$effpop,
      Ne.low = coal_pref$effpop025,
      Ne.high = coal_pref$effpop975))
  
  coal_pref_df$date <- as.Date(convert_decimal_date(tree_mrsd)) - days(round(coal_pref_df$date*365.25,0))
  coal_pref_df = arrange(coal_pref_df, desc(date))
  
  log_breaks <- function(limits) {
    10^pretty(log10(range(limits)))
  }
  
  log_format <- function(x) {
    parse(text = paste("10^", round(log10(x)), sep = ""))
  }
  
  x_min <- as_date(x_limits[1], format = "%Y-%m-%d")
  x_max <- as_date(x_limits[2], format = "%Y-%m-%d")
  
  
  gg_phylo <- ggplot(coal_pref_df, aes(date, Ne)) + 
    geom_ribbon(aes(ymin = Ne.low, ymax = Ne.high), fill = "steelblue", alpha = 0.3) +
    geom_line(col = "black", linewidth=1) +
    scale_x_date(date_breaks = "60 days", date_labels = "%b %Y",
                 limits = c(x_min, x_max)) + 
    scale_y_continuous(trans = "log10", breaks = log_breaks, labels = log_format) +
    ylab("Effective Population Size (Ne)") +
    xlab(" ") +
    labs(fill = "95% CI") + 
    theme_minimal() +
    theme(plot.margin = unit(c(2,0.5,2,0.5),"cm"),
          legend.direction = "vertical",
          legend.position= c(0.9, 0.8), 
          strip.text = element_text(size=26, face="bold"),
          strip.background = element_blank(),
          legend.key.size = unit(2,"line"),
          legend.key.width = unit(1,"line"),
          legend.text = element_text(size=16, face="bold"),
          legend.title = element_text(size=18, face="bold"),
          axis.title.x = element_text(size=24, face="bold"),
          axis.title.y = element_text(size=24, face="bold"),
          axis.text.x = element_text(face="bold", size=12, vjust=1, 
                                     hjust=1, angle=45),
          axis.text.y = element_text(size=12, face="bold"),
          plot.title = element_text(size=28, face="bold"))
  
  return(gg_phylo)
}