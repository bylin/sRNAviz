library(ggplot2)

#his_df = read.table('~/temp/his.tsv')
#colnames(his_df) = c("Position", "Most common base", "Frequency")
#ss = '(((((((,,<<<<________>>>>,<<<<<_______>>>>>,,,,,<<<<<_______>>>>>))))))):'
#bonds=c(rep("astem", 7), rep("dloop", 4), rep("acloop", 5), rep("tloop", 5))

get_coords = function(ss) {
  coords_df = read.table('~/temp/coords.tsv', header=TRUE)
  return(coords_df)
}

merge = function(coords_df, df, outline=NULL, fill_discrete=NULL, fill_continuous=NULL, annotation=NULL) {
  if (!is.null(fill_discrete) & !is.null(fill_continuous)) {
    stop("Fill must either be discrete or continuous, not both")
  }
  if (!is.null(outline)) {
    if (!(outline %in% colnames(df))) {
      stop("outline must be a column in input df")
    }
    coords_df['outline'] = df[outline]
  }
  if (!is.null(fill_discrete)) {
    if (!(fill_discrete %in% colnames(df))) {
      stop("fill_discrete must be a column in input df")
    }
    coords_df['fill'] = df[fill_discrete]
  }
  if (!is.null(fill_continuous)) {
    if (!(fill_continuous %in% colnames(df))) {
      stop("fill_continuous must be a column in input df")
    }
    coords_df['fill'] = df[fill_continuous]
  }
  if (!is.null(annotation)) {
    if (!(annotation %in% colnames(df))) {
      stop("annotation must be a column in input df")
    }
    coords_df['annotation'] = df[annotation]
  }
  return(coords_df)
}

bonds_layer = function(ss, bonds=NULL) {
  base_pairs = list(c(1, 72), c(2, 71), c(3, 70), c(4, 69), c(5, 68), c(6, 67), c(7, 66), c(10, 25), c(11, 24), c(12, 23), c(13, 22), c(27, 43), c(28, 42), c(29, 41), c(30, 40), c(31, 39), c(49, 65), c(50, 64), c(51, 63), c(52, 62), c(53, 61))
  # run positional analysis, get all positions with a ':', grab coords
  ss_df = sapply(base_pairs, function(bp) c(x=coords_df$px[bp[1]], xend=coords_df$px[bp[2]], y=coords_df$py[bp[1]], yend=coords_df$py[bp[2]]))
  ss_df = data.frame(t(ss_df))
  if (!is.null(bonds)) {
    if (length(bonds) != nrow(ss_df)) {
      stop("bonds argument does not match number of bonds")
    }
    ss_df['bonds'] = as.factor(bonds)
    geom = geom_segment(data=ss_df, aes(x=x, y=y, xend=xend, yend=yend, linetype=bonds))
  }
  else {
    geom = geom_segment(data=ss_df, aes(x=x, y=y, xend=xend, yend=yend))
  }
  return(geom)
}

draw_tRNA = function(ss, df, outline=NULL, fill_discrete=NULL, fill_continuous=NULL, annotation=NULL, bonds=NULL, title='', save_file=NULL) {
  # Fold 
  coords_df = get_coords(ss)
  coords_df = merge(coords_df, df, outline=outline, fill_discrete=fill_discrete, fill_continuous=fill_continuous, annotation=annotation)
  
  plot = ggplot() + bonds_layer(ss, bonds)
  
  # outline layer
  if ('outline' %in% colnames(coords_df)) plot = plot + geom_point(data=coords_df, aes(x=px, y=py, color=outline), size=9, shape=1, stroke=1.5) + scale_color_discrete(guide=guide_legend(title=outline))
  else plot = plot + geom_point(data=coords_df, aes(x=px, y=py), color='black', size=9, shape=1, stroke=1.5)
  
  # color fill layer
  if ('fill' %in% colnames(coords_df)) {
    if (!is.null(fill_discrete)) plot = plot + geom_point(data=coords_df, aes(x=px, y=py, fill=fill), shape=21, color='white', size=8.5) + scale_fill_discrete(guide=guide_legend(title=fill_discrete))
    else plot = plot + geom_point(data=coords_df, aes(x=px, y=py, fill=fill), shape=21, color='white', size=8.5) + scale_fill_continuous(guide=guide_legend(title=fill_continuous), low="gray", high="white")
  }
  else plot = plot + geom_point(data=coords_df, aes(x=px, y=py), color='white', size=8.5)
  
  # quality of life improvements
  # remove axes and labels
  plot = plot + theme_bw() + theme(line=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.border=element_blank())
  # maintain aspect ratio for cloverleaf proportions
  #plot = plot + coord_fixed(ratio=1-max(nchar(colnames(coords_df)))/100)
  # title
  plot = plot + ggtitle(title)
  # Save plot to file for proper dimensions
  # ggplot considers the legend as part of the plot, which screws up the cloverleaf proportions. The legend has indeterminate length.
  # This is also not fixable by changing the aspect ratio with coord_fixed, which calculates plot dimensions before adding the legend. Best thing to do is change the output width. 
  # Empirically, a decent estimate is to add (length of longest legend label) / 5. This would add 1.6 inches to the plot width for the string "Most common base".
  width = 7 + max(nchar(colnames(coords_df)))/5
  if (!is.null(save_file)) ggsave(plot, file=save_file, dpi=300, width=width, height=7)
  return(plot)
}

#plot = draw_tRNA(ss, his_df, save_file='plot.png', outline="Most common base", fill_continuous="Frequency", bonds=bonds)
#plot
