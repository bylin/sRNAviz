suppressWarnings(suppressMessages(library(Biostrings)))
library(ggplot2)
library(reshape2)
library(RColorBrewer)
suppressMessages(library(dplyr))
library(stringr)
suppressMessages(library(tidyr))
theme_set(theme_bw())
library(scales)
options(repr.plot.width=7, repr.plot.height=4)
isotypes = c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'iMet', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val')
fills = c('A'='#ffd92f', 'C'='#4daf4a', 'G'='#e41a1c', 'U'='#377eb8', 'A:U'='#93da69', 'U:A'='#93da69', 'G:C'='#c1764a', 'C:G'='#c1764a', 'G:U'='#b26cbd', 'U:G'='#b26cbd', "A:C"='gray30', 'C:A'='gray30', 'A:G'='gray30', 'G:A'='gray30', 'C:U'='gray30', 'U:C'='gray30', 'A:A'='gray30', 'C:C'='gray30', 'G:G'='gray30', 'U:U'='gray30', 'A:-'='gray30', '-:A'='gray30', 'C:-'='gray30', '-:C'='gray30', 'G:-'='gray30', '-:G'='gray30', 'U:-'='gray30', '-:U'='gray30', '-'='gray60', '-:-'='gray60')

load('identities.RData')

calculate_position_specific_scores = function(seq = "", seqname = "", clade = "", isotype = "", anticodon = "") {
  # get subset of tRNAs and write to file
  subset = identities %>% select_('species', 'seqname', 'isotype', 'clade', 'anticodon', 'quality') %>% filter_('quality')
  if (clade != "") subset = subset %>% filter_(paste0("clade == '", clade, "'"))
  if (isotype != "") subset = subset %>% filter_(paste0("isotype == '", isotype, "'"))
  if (isotype != "" & anticodon != "") subset = subset %>% filter_(paste0("anticodon == '", anticodon, "'"))
  if (dim(subset)[1] < 5) return('Rare anticodon; could not generate consensus')
  euk_seqs = readDNAStringSet(filepath = '/projects/lowelab/users/blin/identity/euk-isotypes/euk-tRNAs.fa', format = 'fasta')
  names(euk_seqs) = str_replace(str_extract(names(euk_seqs), '\\S+'), '\\|', '_')
  writeXStringSet(euk_seqs[match(subset$seqname, names(euk_seqs))], filepath = 'subset.fa')

  # create covariance model
  system('cmalign -g --notrunc --matchonly -o subset.sto /projects/lowelab/users/blin/tRNAscan/models/domain-specific/euk-num-092016.cm subset.fa > /dev/null')
  system('cmbuild --hand --enone -F subset.cm subset.sto > /dev/null')

  # remove intron from our tRNA
  #   align our tRNA to the numbering model
  seq = DNAStringSet(seq)
  names(seq) = seqname
  writeXStringSet(seq, filepath = paste0(seqname, "-raw.fa"))
  system(paste('cmalign -g --notrunc --matchonly -o', paste0(seqname, "-raw.sto"), '/projects/lowelab/users/blin/tRNAscan/models/domain-specific/euk-num-092016.cm', paste0(seqname, "-raw.fa"), " > /dev/null"))
  
  #  rewrite tRNA to file and realign to subset model
  seq = str_replace_all(str_extract(as.character(read.delim(paste0(seqname, '-raw.sto'))[2, ]), '[AGCU-]+$'), '-', '')
  seq = RNAStringSet(seq)
  names(seq) = seqname
  writeXStringSet(seq, filepath = paste0(seqname, ".fa"))
  system(paste0('cmalign -g --notrunc --matchonly --tfile ', seqname, '.tfile -o ', seqname, '.sto subset.cm ', seqname, '.fa > /dev/null'))
    
  # parse output
  system(paste0('python parse-parsetree.py ', seqname, '.tfile > ', seqname, '.bits'))
  bits = read.table(paste0(seqname, '.bits'), header = FALSE) %>%
    mutate(Position = factor(V1, c('1:72', '2:71', '3:70', '4:69', '5:68', '6:67', '7:66', '8', '9', '10:25', '11:24', '12:23', '13:22', '14', '15', '16', '17', '18', '19', '20', '20a', '21', '26', '27:43', '28:42', '29:41', '30:40', '31:39', '32', '33', '34', '35', '36', '37', '38', '44', '45', '46', '47', '48', '49:65', '50:64', '51:63', '52:62', '53:61', '54', '55', '56', '57', '58', '59', '60', '73'))) %>%
    mutate(Bits = V2) %>%
    mutate(Identity = V3) %>%
    mutate(Source = "Input") %>%
    select(-V1, -V2) %>%
    filter(!is.na(Position))

  # emit and read in consensus sequence
  # cmemit is not good at putting deletions/insertions into consensus emit (-c switch), so exponentiate instead
  system('cmemit -N 1 --exp 7 -o subset-cons.fa subset.cm > /dev/null')
  system('cmalign -g --notrunc --matchonly --tfile subset-cons.tfile -o subset-cons.sto subset.cm subset-cons.fa > /dev/null')
  system('python parse-parsetree.py subset-cons.tfile > subset-cons.bits')

  bits = rbind(bits, read.table('subset-cons.bits', header = FALSE) %>%
    mutate(Position = factor(V1, c('1:72', '2:71', '3:70', '4:69', '5:68', '6:67', '7:66', '8', '9', '10:25', '11:24', '12:23', '13:22', '14', '15', '16', '17', '18', '19', '20', '20a', '21', '26', '27:43', '28:42', '29:41', '30:40', '31:39', '32', '33', '34', '35', '36', '37', '38', '44', '45', '46', '47', '48', '49:65', '50:64', '51:63', '52:62', '53:61', '54', '55', '56', '57', '58', '59', '60', '73'))) %>%
    mutate(Bits = V2) %>%
    mutate(Identity = V3) %>%
    mutate(Source = "Consensus") %>%
    select(-V1, -V2) %>%
    filter(!is.na(Position)))
    
  # compare consensus and our tRNA; fix instances where cmemit did not output most likely identity
  bits = bits %>% group_by(Position) %>%
    arrange(Source) %>% 
    summarize(Bits = Bits[2] - Bits[1], Consensus = Identity[1], Identity = Identity[2]) %>% 
    ungroup() %>%
    mutate(Consensus = ifelse(Bits > 0, as.character(Identity), as.character(Consensus))) %>%
    mutate(Bits = ifelse(Bits > 0, 0, Bits))

  # clean up
  system('rm subset.sto subset.cm subset.fa')
  system('rm subset-cons.fa subset-cons.bits subset-cons.tfile subset-cons.sto')
  system(paste0('rm ', seqname, '.fa ', seqname, '-raw.fa ', seqname, '.bits ', seqname, '.sto ', seqname, '-raw.sto ', seqname, '.tfile')) 

  return(bits)
}
    
calculate_scores_multiseq = function(seqs = "", clade = "Eukaryota", isotype = "", anticodon = "") {  
  multi_bits = data.frame(Position = character(0), Bits = character(0), Identity = character(0), Consensus = character(0), Clade = character(0), Isotype = character(0), Anticodon = character(0))
  for (i in 1:length(seqs)) {
    bits = calculate_position_specific_scores(seqs[i], names(seqs)[i], clade, isotype, anticodon)
    if (typeof(bits) == "character" && bits == 'Rare anticodon; could not generate consensus') next
    bits = bits %>% mutate(Clade = clade, Isotype = isotype, Anticodon = anticodon, tRNA = names(seqs)[i])
    multi_bits = rbind(multi_bits, bits)
  }
  
  return(multi_bits)
}


trnas = read.table('hg19-tRNAs-unique.tsv', header = FALSE, stringsAsFactors = FALSE)
colnames(trnas) = c('seqname', 'isotype', 'anticodon', 'seq')
isodecoders = trnas %>% select(isotype, anticodon) %>% unique()

for (i in 1:nrow(isodecoders)) {
  current_isotype = isodecoders[i, ]$isotype
  current_anticodon = isodecoders[i, ]$anticodon
  df = trnas %>% filter(isotype == current_isotype, anticodon == current_anticodon)
  df = df %>% mutate(num = as.integer(str_match(seqname, '(\\d+)-1$')[, 2])) %>% arrange(num)
  seqs = df$seq
  names(seqs) = df$seqname
  image_file = paste0('bitcharts/hg19-tRNA-', current_isotype, '-', current_anticodon, '.png')
  print(paste0('Generating ', current_isotype, '-', current_anticodon, ' plot...'))
  if (file.exists(image_file)) {
    print("Image file already exists, skipping")
    next
  }
  bits = calculate_scores_multiseq(seqs, clade = "Mammalia", isotype = current_isotype, anticodon = current_anticodon)
  if (dim(bits)[1] == 0) {
    print("Failed; most likely due to rare anticodon")
    next
  }

  plot = bits %>% select(Position, Bits, Identity, Consensus, tRNA) %>%
    ggplot() + geom_bar(aes(x = Position, y = Bits, fill = Identity), size = 0.1, color = 'gray20', width = 0.8, position = 'dodge', stat = 'identity') + 
    facet_wrap( ~ tRNA, nrow = 1) +
    geom_text(aes(x = Position, color = Identity, label = Identity), size = 3, y = min(bits$Bits) - 0.1 * (max(bits$Bits) - min(bits$Bits))) +
    geom_text(aes(x = Position, color = Consensus, label = Consensus), size = 3, y = 0.1 * (max(bits$Bits) - min(bits$Bits))) +
    coord_flip() +
    scale_fill_manual(values=fills) +
    scale_color_manual(values=fills) + 
    scale_x_discrete(limits = rev(levels(bits$Position))) + 
    scale_y_continuous(limits = c(min(bits$Bits) - 0.15 * (max(bits$Bits) - min(bits$Bits)), 0.15 * (max(bits$Bits) - min(bits$Bits)))) + 
    theme(legend.position='none') + 
    ylab("Score difference between consensus and input (bits)")
  
  ggsave(plot, file = image_file, width = 3 * length(seqs), height = 7, dpi = 250, limitsize = FALSE)
}
