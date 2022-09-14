get_color <- function(levels,this_color_scheme) {
  levels <- full_level_name(levels)
  stopifnot(all(levels %in% scheme_level_names()))
  color_vals <- color_scheme_get(this_color_scheme)[levels]
  unlist(color_vals, use.names = FALSE)
}

full_level_name <- function(x) {
  map <- c(
    l = "light",
    lh = "light_highlight",
    m = "mid",
    mh = "mid_highlight",
    d = "dark",
    dh = "dark_highlight",
    light = "light",
    light_highlight = "light_highlight",
    mid = "mid",
    mid_highlight = "mid_highlight",
    dark = "dark",
    dark_highlight = "dark_highlight"
  )
  unname(map[x])
}

scheme_level_names <- function() {
  c("light",
    "light_highlight",
    "mid",
    "mid_highlight",
    "dark",
    "dark_highlight")
}

mcmc_intervals_compare <- function( model1, model2,
                                    pars = character(),
                                    alpha=0.1, prob=1-alpha/2, prob_outer=1,
                                    model1.name="1", model2.name="2",
                                    useBetas=TRUE,
                                    colorschemes=c("pink","blue"),
                                    colorlist=NULL){
  model1_data <- mcmc_intervals_data(model1, 
                                     pars, prob=prob, prob_outer=prob_outer)
  model1_data$model <- model1.name
  model2_data <- mcmc_intervals_data(model2,
                                     pars, prob=prob, prob_outer=prob_outer)
  model2_data$model <- model2.name
  pars <- unique(model1_data$parameter)
  model1_data$parameter <- paste(model1_data$parameter,"_1",sep="")
  data <- rbind(model1_data,model2_data)
  
  x_lim <- range(c(data$ll, data$hh))
  x_range <- diff(x_lim)
  x_lim[1] <- x_lim[1] - 0.05 * x_range
  x_lim[2] <- x_lim[2] + 0.05 * x_range
  
  # faint vertical line at zero if zero is within x_lim
  layer_vertical_line <- if (0 > x_lim[1] && 0 < x_lim[2]) {
    vline_0(color = "gray90", size = 0.5)
  } else {
    NULL
  }
  
  args_outer <- list(
    mapping = aes_(x = ~ ll, xend = ~ hh, y = ~ parameter, yend = ~ parameter, col = ~ model)#,
    #color = get_color("mid")
  )
  args_inner <- list(
    mapping = aes_(x = ~ l, xend = ~ h, y = ~ parameter, yend = ~ parameter, col = ~ model),
    size = 2,
    show.legend = FALSE
  )
  args_point <- list(
    mapping = aes_(x = ~ m, y = ~ parameter, col = ~ model, fill = ~ model),
    data = data,
    size = 2.5,
    shape = 21
  )
  
  #args_inner$color <- get_color("dark")
  #args_point$color <- get_color("dark_highlight")
  #args_point$fill <- get_color("light")
  
  point_func <- geom_point
  
  layer_outer <- do.call(geom_segment, args_outer)
  layer_inner <- do.call(geom_segment, args_inner)
  layer_point <- do.call(point_func, args_point)
  
  betasTable <- rbind(pars,sapply(1:length(pars),function(i){ bquote(beta[.(i)]) }))
  ordered_pars <- as.character(pars) #sort(as.character(pars))
  if (length(data$parameter) > 2) {
    y_dummies <- paste(ordered_pars[1:(length(pars)-1)],"_",2,sep="")
    y_limits <- rev(c( rbind( paste(ordered_pars,"_",1,sep=""),
                              ordered_pars,
                              c(y_dummies,NA) ) )) #rev(sort(as.character(c(data$parameter,y_dummies))))
    y_limits <- y_limits[!is.na(y_limits)]
    if (useBetas){     
      y_labels_pars <- betasTable[2,]
    } else {
      y_labels_pars <- ordered_pars }
    y_labels <- rev(c(rbind(y_labels_pars,"",""))[1:length(y_limits)])
  } else {
    y_dummies <- data$parameter[length(data$parameter)]
    y_limits <- rev(sort(as.character(c(data$parameter))))
    y_labels <- c("",y_dummies)
  }
  
  if (is.null(colorlist)){
    colorlist[1] <- get_color("mid",colorschemes[1]) # mid1
    colorlist[2] <- get_color("mid",colorschemes[2]) # mid2
    colorlist[3] <- get_color("light",colorschemes[1]) # lt1
    colorlist[4] <- get_color("light",colorschemes[2]) # lt2
  }
  
  ggplot(data) +
    layer_vertical_line +
    layer_outer +
    layer_inner +
    layer_point +
    scale_color_manual(values=colorlist[1:2]) +
    scale_fill_manual(values=colorlist[3:4]) +
    #scale_color +
    #scale_fill +
    scale_y_discrete(limits = y_limits,
                     labels = y_labels ) +
    xlim(x_lim) +
    bayesplot_theme_get() +
    legend_move("bottom") +
    yaxis_text(face = "bold") +
    yaxis_title(FALSE) +
    #yaxis_ticks(size = 1) +
    theme(axis.ticks.y = element_blank()) +
    xaxis_title(FALSE)
}