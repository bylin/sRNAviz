#!/usr/bin/env Rscript
library(ggplot2)
library(RColorBrewer)
suppressMessages(library(dplyr))
library(stringr)
theme_set(theme_bw())

# Setup convenience globals
isotypes = c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'iMet', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val')

paired_positions = c('X1.72'='(A) 1:72', 'X2.71'='(A) 2:71', 'X3.70'='(A) 3:70', 'X4.69'='(A) 4:69', 'X5.68'='(A) 5:68', 'X6.67'='(A) 6:67', 'X7.66'='(A) 7:66', 'X8.14'='(A) 8:14', 'X9.23'='(3) 9:23', 'X10.25'='(D) 10:25', 'X10.45'='(3) 10:45', 'X11.24'='(D) 11:24', 'X12.23'='(D) 12:23', 'X13.22'='(D) 13:22', 'X15.48'='(3) 15:48', 'X18.55'='(3) 18:55', 'X19.56'='(3) 19:56', 'X22.46'='(3) 22:46', 'X26.44'='(C) 26:44', 'X27.43'='(C) 27:43', 'X28.42'='(C) 28:42', 'X29.41'='(C) 29:41', 'X30.40'='(C) 30:40', 'X31.39'='(C) 31:39', 'X49.65'='(T) 49:65', 'X50.64'='(T) 50:64', 'X51.63'='(T) 51:63', 'X52.62'='(T) 52:62', 'X53.61'='(T) 53:61', 'X54.58'='(3) 54:58')
paired_identities = c('GC', 'AU', 'UA', 'CG', 'GU', 'UG', 'PairDeletion', 'PurinePyrimidine', 'PyrimidinePurine', 'StrongPair', 'WeakPair', 'Wobble', 'Paired', 'Bulge', 'Mismatched')
paired_colors = c('GC'='gray20', 'AU'='gray20', 'UA'='gray20', 'CG'='gray20', 'GU'='gray20', 'UG'='gray20', 'PairDeletion'='gray20', 'PurinePyrimidine'='gray40', 'PyrimidinePurine'='gray40', 'StrongPair'='gray40', 'WeakPair'='gray40', 'Wobble'='gray40', 'Paired'='gray40', 'Bulge'='gray40', 'Mismatched'='gray40')

single_positions = c('X8'='8', 'X9'='9', 'X14'='(D) 14', 'X15'='(D) 15', 'X16'='(D) 16', 'X17'='(D) 17', 'X17a'='(D) 17a', 'X18'='(D) 18', 'X19'='(D) 19', 'X20'='(D) 20', 'X20a'='(D) 20a', 'X20b'='(D) 20b', 'X21'='(D) 21', 'X26'='26', 'X32'='(C) 32', 'X33'='(C) 33', 'X34'='(C) 34', 'X35'='(C) 35', 'X36'='(C) 36', 'X37'='(C) 37', 'X38'='(C) 38', 'X44'='(V) 44', 'X45'='(V) 45', 'X46'='(V) 46', 'X47'='(V) 47', 'X48'='(V) 48', 'X54'='(T) 54', 'X55'='(T) 55', 'X56'='(T) 56', 'X57'='(T) 57', 'X58'='(T) 58', 'X59'='(T) 59', 'X60'='(T) 60', 'X73'='73')
single_identities = c('A', 'C', 'G', 'U', 'Deletion', 'Purine', 'Pyrimidine', 'Weak', 'Strong', 'Amino', 'Keto', 'B', 'D', 'H', 'V')
single_colors = c('A'='gray20', 'C'='gray20', 'G'='gray20', 'U'='gray20', 'Deletion'='gray20', 'Purine'='gray40', 'Pyrimidine'='gray40', 'Weak'='gray40', 'Strong'='gray40', 'Amino'='gray40', 'Keto'='gray40', 'B'='gray40', 'D'='gray40', 'H'='gray40', 'V'='gray40')

# Get user input
args = commandArgs(trailingOnly=TRUE)
seq = args[1]
input_species_clade = args[2]
input_isotypes = args[3]
session_id = args[4]
input_cutoff = args[5]
input_isotypes = unlist(str_split(input_isotypes, ','))
single_plot_path = paste0('/tmp/single-plot-', session_id, '.png')
paired_plot_path = paste0('/tmp/paired-plot-', session_id, '.png')

# Read in data
freqs = read.table('single-trna-heatmap/freqs.tsv', sep='\t', header=TRUE, stringsAsFactors=FALSE)
freqs = freqs %>% filter(cutoff == as.numeric(input_cutoff)) %>% select(-cutoff)
genome_table = read.delim('single-trna-heatmap/genome_table+.txt', header=FALSE, sep='\t', stringsAsFactors=FALSE)
if (input_species_clade %in% genome_table$V5) {
  input_clade = input_species_clade
} else {
  input_clade = genome_table[genome_table$V2 == input_species_clade, ]$V5
}

# write fasta file
fasta_file = paste0('/tmp/seq-', session_id, '.fa')
fasta_handle = file(fasta_file) 
writeLines(c(paste0(">", fasta_file), seq), fasta_file)
close(fasta_handle)

# Align seq
alignment_file = paste0('/tmp/alignment-', session_id, '.sto')
model = '/projects/lowelab/users/blin/tRNAscan/models/domain-specific/euk-num-092016.cm'
output = system(paste('cmalign -g --notrunc', model, fasta_file, '| tee', alignment_file), intern=TRUE)
seq = tail(unlist(str_split(output[4], '\\s+')), 1)
ss = tail(unlist(str_split(output[6], '\\s+')), 1)

# Get seq numbering
output = system(paste0('printf "', seq, '\\n', ss, '" | python single-trna-heatmap/position_interface.py'), intern=TRUE)
df = data.frame(t(matrix(unlist(str_split(output, '\\s')), nrow=2)), "Input", "Input", stringsAsFactors=FALSE)
codes = c("A"="A", "C"="C", "G"="G", "U"="U", "-"="Deletion", "."="Deletion", "A:A"="Mismatched", "G:G"="Mismatched", "C:C"="Mismatched", "U:U"="Mismatched", "A:G"="Mismatched", "A:C"="Mismatched", "C:A"="Mismatched", "C:U"="Mismatched", "G:A"="Mismatched", "U:C"="Mismatched", "A:-"="Bulge", "U:-"="Bulge", "C:-"="Bulge", "G:-"="Bulge", "-:A"="Bulge", "-:G"="Bulge", "-:C"="Bulge", "-:U"="Bulge", "A:U"="AU", "U:A"="UA", "C:G"="CG", "G:C"="GC", "G:U"="GU", "U:G"="UG", "-:-"="PairDeletion")
df$identity = codes[df$X2]

# Filter out varm and insertions
df = df %>% filter(!str_detect(X1, "i") & !str_detect(X1, 'V'))

# Wrangle to fit main data frame
df$X1 = gsub("^", "X", gsub(':', '.', df$X1))
df = df[, c(3, 1, 5, 4)]
colnames(df) = c("isotype", "positions", "identity", "clade")
freqs = rbind(freqs, df)

# Wrange main data frame for plotting
matches_input = function(isotype, positions, identity) {
  codes = list(A='A', C='C', G='G', U='U', Deletion='Deletion', Purine=c('A', 'G', 'Purine'), Pyrimidine=c('C', 'U', 'Pyrimidine'), Weak=c('A', 'U', 'Weak'), Strong=c('G', 'C', 'Strong'), Amino=c('A', 'C', 'Amino'), Keto=c('G', 'U', 'Keto'), B=c('C', 'G', 'U', 'B', 'Strong', 'Pyrimidine', 'Keto'), D=c('A', 'G', 'U', 'D', 'Purine', 'Weak', 'Keto'), H=c('A', 'C', 'U', 'H', 'Amino', 'Weak', 'Pyrimidine'), V=c('A', 'C', 'G', 'V', 'Amino', 'Purine', 'Strong'), GC='GC', AU='AU', UA='UA', CG='CG', GU='GU', UG='UG', PairDeletion='PairDeletion', PurinePyrimidine=c('AU', 'GC', 'PurinePyrimidine'), PyrimidinePurine=c('UA', 'CG', 'PyrimidinePurine'), StrongPair=c('GC', 'CG', 'StrongPair'), WeakPair=c('AU', 'UA', 'WeakPair'), Wobble=c('GU', 'UG', 'Wobble'), Paired=c('AU', 'UA', 'CG', 'GC', 'GU', 'UG', 'Paired', 'PurinePyrimidine', 'PyrimidinePurine', 'StrongPair', 'WeakPair', 'Wobble'), Bulge=c('Bulge'), Mismatched=c('AA', 'GG', 'CC', 'UU', 'AG', 'AC', 'CA', 'CU', 'GA', 'UC', 'Mismatched', 'Paired', 'PurinePyrimidine', 'PyrimidinePurine', 'StrongPair', 'WeakPair', 'Wobble'))
  seq_identity = freqs[freqs$positions == positions & freqs$clade == "Input", ]$identity
  if (length(seq_identity) != 1) print(paste("Too many codes", isotype, positions, identity))
  return(ifelse(seq_identity == as.character(identity),
                "Match", 
                ifelse(seq_identity %in% codes[[as.character(identity)]], "Subset", "Conflict")))
}
freqs = freqs %>% 
  filter(isotype %in% c("Input", input_isotypes) & clade %in% c("Eukarya", "Input", input_clade)) %>% 
  mutate(category=ifelse(clade != "Input", paste0(isotype, ' - ', clade), "Your Sequence")) %>%
  rowwise() %>% mutate(match=matches_input(isotype, positions, identity))

# Write data to file
write.table(freqs[, c('isotype', 'positions', 'identity', 'clade', 'match')], file=paste0('/tmp/identities-', session_id, '.tsv'), quote=FALSE, sep='\t', row.names=FALSE)

# Single plot
plot = freqs %>%
  filter(positions %in% names(single_positions)) %>%
  mutate(positions=factor(positions, names(single_positions))) %>%
  mutate(identity=factor(identity, single_identities)) %>%
  ggplot() + geom_tile(aes(x=positions, y=category, fill=identity, color=identity), width=0.9, height=0.9, size=0.5) + 
    geom_point(aes(x=positions, y=category, shape=match)) +
    scale_shape_manual(values=c(4, 1, 2), labels=c("Conflict", "Match", "Subset")) +
    scale_x_discrete(labels=single_positions) +
    scale_color_manual(values=single_colors) +
    scale_fill_manual(values=c(brewer.pal(5, "Set1"), brewer.pal(12, "Set3"))) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position='bottom', legend.text=element_text(size=8)) + 
    guides(fill=guide_legend(title=NULL, nrow=2, keywidth=0.8, keyheight=0.8, order=2), 
           color=guide_legend(title=NULL, nrow=2, keywidth=0.8, keyheight=0.8, order=2), 
           shape=guide_legend(title=NULL, keywidth=0.8, keyheight=0.8, order=1)) + 
    xlab('Position') + ylab('Dataset')
ggsave(plot, file=single_plot_path, width=8, height=1.65+0.42*length(input_isotypes))

# Paired plot
plot = freqs %>% 
  filter(positions %in% names(paired_positions)) %>%
  mutate(positions=factor(positions, names(paired_positions))) %>%
  mutate(identity=factor(identity, paired_identities)) %>%
  ggplot() + geom_tile(aes(x=positions, y=category, fill=identity, color=identity), width=0.9, height=0.9, size=0.5) + 
    geom_point(aes(x=positions, y=category, shape=match)) +
    scale_shape_manual(values=c(4, 1, 2), labels=c("Conflict", "Match", "Subset")) +
    scale_x_discrete(labels=paired_positions) +
    scale_color_manual(values=paired_colors) +
    scale_fill_manual(values=c(brewer.pal(5, "Set1"), brewer.pal(12, "Set3"))) +
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.position='bottom', legend.text=element_text(size=8)) + 
    guides(fill=guide_legend(title=NULL, nrow=2, keywidth=0.8, keyheight=0.8, order=2),
           color=guide_legend(title=NULL, nrow=2, keywidth=0.8, keyheight=0.8, order=2),
           shape=guide_legend(title=NULL, keywidth=0.8, keyheight=0.8, order=1)) + 
    xlab('Position') + ylab('Dataset')
ggsave(plot, file=paired_plot_path, width=8, height=1.75+0.48*length(input_isotypes))

# clean up fasta file
system(paste0("rm ", fasta_file))