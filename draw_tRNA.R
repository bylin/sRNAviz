library(ggplot2)

his_df = read.table('~/temp/his.tsv')
colnames(his_df) = c("Position", "Most common base", "Frequency")
ss = '(((((((,,<<<<________>>>>,<<<<<_______>>>>>,,,,,<<<<<_______>>>>>))))))):'
bonds=c(rep("astem", 7), rep("dloop", 4), rep("acloop", 5), rep("tloop", 5))

fold = function(seq) return(seq)

get_coords = function(ss) {
  coords_df = read.table('~/temp/coords.tsv', header=TRUE)
  return(coords_df)
}

merge = function(coords_df, df, outline=NULL, fill_discrete=NULL, fill_continuous=NULL, annotation=NULL) {
  if (nrow(df) != 73) stop("Input df must have the same number of rows as the length of the input sequence")
  if (!is.null(fill_discrete) & !is.null(fill_continuous)) stop("Fill must either be discrete or continuous, not both")
  if (!is.null(outline)) {
    if (!(outline %in% colnames(df))) stop(paste(outline, "must be a column in input df"))
    coords_df['outline'] = df[outline]
  }
  if (!is.null(fill_discrete)) {
    if (!(fill_discrete %in% colnames(df))) stop(paste(fill_discrete, "must be a column in input df"))
    coords_df['fill'] = df[fill_discrete]
  }
  if (!is.null(fill_continuous)) {
    if (!(fill_continuous %in% colnames(df))) stop(paste(fill_continuous, "must be a column in input df"))
    coords_df['fill'] = df[fill_continuous]
  }
  if (!is.null(annotation)) {
    if (!(annotation %in% colnames(df))) stop(paste(annotation, "must be a column in input df"))
    if (typeof(unlist(df[annotation])) != "character") stop("annotation must be of type string")
    coords_df['annotation'] = df[annotation]
  }
  return(coords_df)
}

add_bonds_layer = function(plot, ss, coords_df, bonds=NULL) {
  base_pairs = list(c(1, 72), c(2, 71), c(3, 70), c(4, 69), c(5, 68), c(6, 67), c(7, 66), c(10, 25), c(11, 24), c(12, 23), c(13, 22), c(27, 43), c(28, 42), c(29, 41), c(30, 40), c(31, 39), c(49, 65), c(50, 64), c(51, 63), c(52, 62), c(53, 61))
  # run positional analysis, get all positions with a ':', grab coords
  ss_df = sapply(base_pairs, function(bp) c(x=coords_df$px[bp[1]], xend=coords_df$px[bp[2]], y=coords_df$py[bp[1]], yend=coords_df$py[bp[2]]))
  ss_df = data.frame(t(ss_df))
  if (!is.null(bonds)) {
    if (length(bonds) != nrow(ss_df)) {
      stop("bonds argument does not match number of bonds")
    }
    ss_df['bonds'] = as.factor(bonds)
    plot = plot + geom_segment(data=ss_df, aes(x=x, y=y, xend=xend, yend=yend, linetype=bonds))
  }
  else {
    plot = plot + geom_segment(data=ss_df, aes(x=x, y=y, xend=xend, yend=yend))
  }
  return(plot)
}

add_outline_layer = function(plot, coords_df, outline) {
  grain = 100
  angle = seq(-pi, pi, length = grain)
  df = data.frame()
  if ('outline' %in% colnames(coords_df)) {
    for (i in 1:nrow(coords_df)) {
      row = coords_df[i, ]
      df = rbind(df, data.frame(position=rep(row$position, grain), x=row$px - sin(angle)/50, y=row$py - cos(angle)/50, color=rep(row$outline, grain)))
    }
    df$color = as.factor(df$color)
    plot = plot + geom_path(data=df, aes(x=x, y=y, group=position, color=color), size=1) + scale_color_discrete(guide=guide_legend(title=outline))
  } else if (all(c('A', 'C', 'G', 'T') %in% colnames(coords_df))) {
    for (i in 1:nrow(coords_df)) {
      row = coords_df[i, ]
      df = rbind(df, data.frame(position=rep(row$position, grain), x=row$px - sin(angle)/50, y=row$py - cos(angle)/50, color=ifelse(seq(1,100)/100<row$A, 'A', 'T')))
    }
    plot = plot +geom_path(data=df, aes(x=x, y=y, group=position, color=color), size=1) + scale_color_discrete(guide=guide_legend(title=outline))
  } else {
    for (i in 1:nrow(coords_df)) {
      row = coords_df[i, ]
      df = rbind(df, data.frame(position=rep(row$position, grain), x=row$px - sin(angle)/50, y=row$py - cos(angle)/50))
    }
    plot = plot +geom_path(data=df, aes(x=x, y=y, group=position), color='black', size=1)
  }
  
  return(plot)
}

draw_tRNA = function(seq, df, outline=NULL, fill_discrete=NULL, fill_continuous=NULL, annotation=NULL, bonds=NULL, title='', save_file=NULL) {
  # Fold and map folded secondary structure to x and y coordinates
  ss = fold(seq)
  coords_df = get_coords(ss)
  coords_df = merge(coords_df, df, outline=outline, fill_discrete=fill_discrete, fill_continuous=fill_continuous, annotation=annotation)
  
  plot = ggplot()
  plot = add_bonds_layer(plot, ss, coords_df, bonds)
  plot = add_outline_layer(plot, coords_df, outline)
  
  # color fill layer
  if ('fill' %in% colnames(coords_df)) {
    if (!is.null(fill_discrete)) plot = plot + geom_point(data=coords_df, aes(x=px, y=py, fill=fill), shape=21, color='white', size=8) + scale_fill_discrete(guide=guide_legend(title=fill_discrete))
    else plot = plot + geom_point(data=coords_df, aes(x=px, y=py, fill=fill), color='white', shape=21, size=8) + scale_fill_gradient(guide=guide_legend(title=fill_continuous), low="darkgray", high="white", limits=c(0,1))
  }
  else plot = plot + geom_point(data=coords_df, aes(x=px, y=py), color='white', size=8)
  
  # text annotation layer
  if ('annotation' %in% colnames(coords_df)) {
    size = 10/mean(nchar(coords_df$annotation))
    plot = plot + geom_text(data=coords_df, aes(x=px, y=py, label=annotation), size=size)
  }
  else plot = plot + geom_text(data=coords_df, aes(x=px, y=py, label=position))
  
  # quality of life improvements
  plot = plot + ggtitle(title)
  # remove axes and labels
  plot = plot + theme_bw() + theme(line=element_blank(), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), panel.border=element_blank())
  # Save plot to file for proper dimensions
  # ggplot considers the legend as part of the plot, which can screw up the cloverleaf proportions. This is not fixable by changing the aspect ratio with coord_fixed, which calculates plot dimensions before adding the legend. Best thing to do is change the output width. 
  # Empirically, a decent estimate is to add (length of longest legend label) / 5. This would add 1.6 inches to the plot width for the string "Most common base".
  width = 8 + max(nchar(colnames(coords_df)))/10
  if (!is.null(save_file)) ggsave(plot, file=save_file, dpi=300, width=width, height=8)
  return(plot)
}

#his_df$plz=paste0(his_df$Position, sapply(as.integer(runif(73) * 5), function(x) paste(rep('a', x), collapse='')))
#
#plot = draw_tRNA(ss, his_df, save_file='plot.png', outline="Most common base", fill_continuous="Frequency", bonds=bonds, annotation='plz')
#plot
#
#plot3 = ggplot() + geom_path(data=df, aes(x=x, y=y, group=position, color=color), size=1)
#plot3 = plot3 + geom_point(data=coords_df, aes(x=px, y=py), color='white', shape=21, size=7)
#plot3 = plot3 + geom_text(data=coords_df, aes(x=px, y=py, label=position))
#ggsave(plot3, file='~/temp/plot.png', width=8, height=8, dpi=300)
#
#outline_layer(coords_df, outline=NULL)

