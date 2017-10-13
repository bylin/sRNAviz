suppressWarnings(suppressMessages(library(Biostrings)))
suppressMessages(library(tidyverse))
library(RColorBrewer)
library(stringr)
theme_set(theme_bw())
suppressMessages(library(scales))
options(repr.plot.width=7, repr.plot.height=4)
isotypes = c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'iMet', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val')
positions = c('X1.72'='1:72', 'X2.71'='2:71', 'X3.70'='3:70', 'X4.69'='4:69', 'X5.68'='5:68', 'X6.67'='6:67', 'X7.66'='7:66', 'X8'='8', 'X9'='9', 'X10.25'='10:25', 'X11.24'='11:24', 'X12.23'='12:23', 'X13.22'='13:22', 'X14'='14', 'X15'='15', 'X16'='16', 'X17'='17', 'X17a'='17a', 'X18'='18', 'X19'='19', 'X20'='20', 'X20a'='20a', 'X20b'='20b', 'X21'='21', 'X26'='26','X27.43'='27:43', 'X28.42'='28:42', 'X29.41'='29:41', 'X30.40'='30:40', 'X31.39'='31:39', 'X32'='32', 'X33'='33', 'X34'='34', 'X35'='35', 'X36'='36', 'X37'='37', 'X38'='38', 'X44'='44', 'X45'='45', 'X46'='46', 'X47'='47', 'X48'='48', 'X49.65'='49:65', 'X50.64'='50:64', 'X51.63'='51:63', 'X52.62'='52:62', 'X53.61'='53:61', 'X54'='54', 'X55'='55', 'X56'='56', 'X57'='57', 'X57'='58', 'X59'='59', 'X60'='60', 'X73'='73')
fills = c('A'='#ffd92f', 'C'='#4daf4a', 'G'='#e41a1c', 'U'='#377eb8', 'A:U'='#93da69', 'U:A'='#93da69', 'G:C'='#c1764a', 'C:G'='#c1764a', 'G:U'='#b26cbd', 'U:G'='#b26cbd', "A:C"='gray30', 'C:A'='gray30', 'A:G'='gray30', 'G:A'='gray30', 'C:U'='gray30', 'U:C'='gray30', 'A:A'='gray30', 'C:C'='gray30', 'G:G'='gray30', 'U:U'='gray30', 'A:-'='gray30', '-:A'='gray30', 'C:-'='gray30', '-:C'='gray30', 'G:-'='gray30', '-:G'='gray30', 'U:-'='gray30', '-:U'='gray30', '-'='gray60', '-:-'='gray60', 'N'='gray60', 'N:N'='gray60')

load('identities.RData')
load('clade-isotype-specific.RData')
load('isotype-specific.RData')
load('consensus-IDEs.RData')

calculate_position_specific_scores = function(seq = "", seqname = "", clade = "", isotype = "", anticodon = "") {
  # get subset of tRNAs and write to file
  subset = identities %>% select_('species', 'seqname', 'isotype', 'clade', 'anticodon', 'quality') %>% filter_('quality')
  if (clade != "") subset = subset %>% filter_(paste0("clade == '", clade, "'"))
  if (isotype != "") subset = subset %>% filter_(paste0("isotype == '", isotype, "'"))
  if (anticodon != "") subset = subset %>% filter_(paste0("anticodon == '", anticodon, "'"))
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
    mutate(Position = factor(V1, c('1:72', '2:71', '3:70', '4:69', '5:68', '6:67', '7:66', '8', '9', '10:25', '11:24', '12:23', '13:22', '14', '15', '16', '17', '17a', '18', '19', '20', '20a', '20b', '21', '26', '27:43', '28:42', '29:41', '30:40', '31:39', '32', '33', '34', '35', '36', '37', '38', '44', '45', '46', '47', '48', '49:65', '50:64', '51:63', '52:62', '53:61', '54', '55', '56', '57', '58', '59', '60', '73'))) %>%
    mutate(Bits = V2) %>%
    mutate(Identity = V3) %>%
    mutate(Source = "Input") %>%
    select(-V1, -V2, -V3) %>%
    filter(!is.na(Position))

  # generate consensus sequence. This is a multistep process.
  #   Create a new model that doesn't conform to numbering model. This is important because cmemit can't emit a consensus with gaps.
  system('cmemit --exp 5 -N 1000 -a subset.cm > subset-free.sto')
  system('cmbuild --enone -F subset-free.cm subset-free.sto > /dev/null')
  #   Next, emit and read in consensus sequence
  #   replace lowercase residues with uppercase. Used to replace with Ns, which would lead to the consensus scoring lower than 0. 
  system('cmemit -c subset-free.cm | perl -npe "if(/^[acguACGU]/){s/uc($1)/[agcu]/g}" > subset-cons.fa')
  system('cmalign -g --notrunc --matchonly --tfile subset-cons.tfile -o subset-cons.sto subset.cm subset-cons.fa > /dev/null')
  system('python parse-parsetree.py subset-cons.tfile > subset-cons.bits')

  bits = rbind(bits, read.table('subset-cons.bits', header = FALSE) %>%
    mutate(Position = factor(V1, c('1:72', '2:71', '3:70', '4:69', '5:68', '6:67', '7:66', '8', '9', '10:25', '11:24', '12:23', '13:22', '14', '15', '16', '17', '17a', '18', '19', '20', '20a', '20b', '21', '26', '27:43', '28:42', '29:41', '30:40', '31:39', '32', '33', '34', '35', '36', '37', '38', '44', '45', '46', '47', '48', '49:65', '50:64', '51:63', '52:62', '53:61', '54', '55', '56', '57', '58', '59', '60', '73'))) %>%
    mutate(Bits = V2) %>%
    mutate(Identity = V3) %>%
    mutate(Source = "Consensus") %>%
    select(-V1, -V2, -V3) %>%
    filter(!is.na(Position)))

  # compare consensus and our tRNA
  bits = bits %>% group_by(Position) %>%
    arrange(Source) %>% 
    summarize(Bits = Bits[2] - Bits[1], Consensus = Identity[1], Identity = Identity[2]) %>% 
    ungroup()    

  # clean up
  system('rm subset.sto subset.cm subset.fa')
  system('rm subset-cons.fa subset-cons.bits subset-cons.tfile subset-cons.sto')
  system('rm subset-free.sto subset-free.cm')
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

# helper function for extracting consensus columns
consensus_cols = function(bits) {
  bold_identities = c("A"="A", "C"="C", "G"="G", "U"="U", "AU"="A:U", "UA"="U:A", "GC"="G:C", "CG"="C:G", "GU"="G:U", "UG"="U:G", "Purine"="R", "Pyrimidine"="Y", "PurinePyrimidine"="R:Y", "PyrimidinePurine"="Y:R", "Mismatched"="N N", "Absent"="-")
  isotype = (bits %>% select(Isotype) %>% arrange(desc(Isotype)) %>% unlist %>% unname)[1]
  clade = (bits %>% select(Clade) %>% arrange(desc(Clade)) %>% unlist %>% unname)[1]

  if (isotype == "") df = consensus
  else if (clade %in% c("", "Eukaryota")) df = isotype_specific %>% filter_(paste0("isotype == '", isotype, "'"))
  else df = clade_isotype_specific %>% filter_(paste0("isotype == '", isotype, "'"), paste0("clade == '", clade, "'"))

  features = data.frame(Position = character(0), Identity = character(0), Color = character(0))
  for (position in sapply(str_extract_all(string = bits$Position %>% unique, pattern = "\\d+[ab]?"), function(x) paste0('X', paste(x, collapse = '.')))) {
    identity = df %>% filter(positions == position)
    if (nrow(identity) == 0) features = rbind(features, data.frame(Position = position, Identity = '', Color = "gray30"))
    else {
      identity = identity %>% ungroup %>% select(identity) %>% unlist %>% unname
      if (identity %in% names(bold_identities)) features = rbind(features, data.frame(Position = position, Identity = bold_identities[identity], Color = 'steelblue'))
      else features = rbind(features, data.frame(Position = position, Identity = '', Color = "gray30"))
    }
  }

  features %>% mutate(Position = str_replace(str_replace(Position, "\\.", ":"), "X", "")) %>%
      mutate(Color = as.character(Color)) %>%
      mutate(Bits = 0, Consensus = '', Clade = '', Anticodon = '', Isotype = '', tRNA = 'Consensus')
}

plot_bitchart = function(bits) {
  features = consensus_cols(bits)
  bits %>% mutate(Bits = ifelse(Bits > 0, 0, ifelse(Bits < -15, -15, Bits))) %>% # normalize tRNA scores to 0 (max) and -15 (min)
    rbind(features %>% select(-Color)) %>%
    ggplot() + geom_tile(aes(x = Position, y = tRNA, fill = Bits, alpha = -Bits)) +
      geom_text(aes(x = Position, y = tRNA, label = Identity, color = ifelse(Bits < -5, "white", "black")), size = 3) +
      scale_fill_gradientn(colors = c("mediumpurple4", "firebrick", "white"),
                           values = c(0, 2/3, 1),
                           limits = c(-15, 0)) +
      scale_color_manual(values = c("black" = "black", "white" = "white")) +
      scale_alpha(range = c(0.4, 1)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = features$Color), legend.position = "bottom") + 
      guides(alpha = FALSE, color = FALSE, fill = guide_colorbar(title = "Score", barwidth = 10)) + 
      labs(y = "") + 
      coord_equal()
}

trnas = read.table('hg19-tRNAs-unique.tsv', header = FALSE, stringsAsFactors = FALSE)
colnames(trnas) = c('seqname', 'isotype', 'anticodon', 'seq')
isodecoders = trnas %>% select(isotype, anticodon) %>% unique()

for (i in 1:nrow(isodecoders)) {
  current_isotype = isodecoders[i, ]$isotype
  current_anticodon = isodecoders[i, ]$anticodon
  if (current_isotype != "Undet") {
    df = trnas %>% filter(isotype == current_isotype, anticodon == current_anticodon)
    df = df %>% mutate(num = as.integer(str_match(seqname, '(\\d+)-1$')[, 2])) %>% arrange(num)
  }
  else {
    df = trnas %>% filter(isotype %in% c("Undet", "Sup") | str_detect(seqname, "Und"))
  }
  seqs = df$seq
  names(seqs) = df$seqname
  image_file = paste0('bitcharts/hg19-tRNA-', current_isotype, '-', current_anticodon, '.png')
  print(paste0('Generating ', current_isotype, '-', current_anticodon, ' plot...'))
  if (file.exists(image_file)) {
    print("Image file already exists, skipping")
    next
  }

  if (current_isotype != "Undet") bits = calculate_scores_multiseq(seqs, clade = "Mammalia", isotype = current_isotype, anticodon = current_anticodon)
  else bits = calculate_scores_multiseq(seqs, clade = "Mammalia", anticodon = current_anticodon)
  
  if (dim(bits)[1] == 0) {
    print("Failed; most likely due to rare anticodon")
    next
  }

  plot = bits %>% rowwise() %>% mutate(Identity = ifelse(Consensus == Identity, "", paste(Identity)), Bits = ifelse(Identity == "", 0, Bits)) %>%
    rbind(bits %>% 
            filter(tRNA == tRNA[1]) %>%
            mutate(Identity = Consensus, Bits = 0, tRNA = "Modal feature")) %>%
    mutate(tRNA = factor(tRNA, levels = rev(c("Modal feature", unique(bits$tRNA))))) %>%
    plot_bitchart + theme(axis.text.y = element_text(color = c("red", "gray30", "gray30", "gray30", "gray30")))
  
  ggsave(plot, file = image_file, width = 16, height = 1.7 + 0.26 * (length(seqs) + 2), limitsize = FALSE)
}
