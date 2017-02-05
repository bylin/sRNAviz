library(reshape2)
library(dplyr)
library(stringr)
library(tidyr)
isotypes = c('Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Gln', 'Glu', 'Gly', 'His', 'Ile', 'iMet', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val')

resolve_code = function(codes) {
  codes = unique(codes)
  x = c(A = all(codes %in% 'A'),
        C = all(codes %in% 'C'),
        G = all(codes %in% 'G'),
        U = all(codes %in% 'U'),
        Deletion = all(codes %in% 'Deletion'),
        Purine = all(codes %in% c('A', 'G', 'Purine')),
        Pyrimidine = all(codes %in% c('C', 'U', 'Pyrimidine')),
        Weak = all(codes %in% c('A', 'U', 'Weak')),
        Strong = all(codes %in% c('G', 'C', 'Strong')),
        Amino = all(codes %in% c('A', 'C', 'Amino')),
        Keto = all(codes %in% c('G', 'U', 'Keto')),
        B = all(codes %in% c('C', 'G', 'U', 'B', 'Strong', 'Pyrimidine', 'Keto')),
        D = all(codes %in% c('A', 'G', 'U', 'D', 'Purine', 'Weak', 'Keto')),
        H = all(codes %in% c('A', 'C', 'U', 'H', 'Amino', 'Weak', 'Pyrimidine')),
        V = all(codes %in% c('A', 'C', 'G', 'V', 'Amino', 'Purine', 'Strong')),
        GC = all(codes %in% 'GC'),
        AU = all(codes %in% 'AU'),
        UA = all(codes %in% 'UA'),
        CG = all(codes %in% 'CG'),
        GU = all(codes %in% 'GU'),
        UG = all(codes %in% 'UG'),
        PairDeletion = all(codes %in% 'PairDeletion'), 
        PurinePyrimidine = all(codes %in% c('AU', 'GC', 'PurinePyrimidine')),
        PyrimidinePurine = all(codes %in% c('UA', 'CG', 'PyrimidinePurine')),
        StrongPair = all(codes %in% c('GC', 'CG', 'StrongPair')),
        WeakPair = all(codes %in% c('AU', 'UA', 'WeakPair')),
        Wobble = all(codes %in% c('GU', 'UG', 'Wobble')),
        Paired = all(codes %in% c('AU', 'UA', 'CG', 'GC', 'GU', 'UG', 'Paired', 'PurinePyrimidine', 'PyrimidinePurine', 'StrongPair', 'WeakPair', 'Wobble')),
        Bulge = all(codes %in% 'Bulge'),
        Mismatched = all(codes %in% c('AA', 'GG', 'CC', 'UU', 'AG', 'AC', 'CA', 'CU', 'GA', 'UC', 'Mismatched')))
  return(names(x[which(x)]))
}

print("Reading in identities.tsv")
identities = read.delim('identities.tsv', sep='\t', stringsAsFactors=FALSE)
identities$quality = as.logical(identities$quality)
identities$restrict = as.logical(identities$restrict)
positions = colnames(identities)[which(str_detect(colnames(identities), "X\\d+\\.\\d+$"))]
positions = c(positions, 'X8', 'X9', 'X14', 'X15', 'X16', 'X17', 'X18', 'X19', 'X20', 'X20a', 'X21', 'X26', 'X32', 'X33', 'X34', 'X35', 'X36', 'X37', 'X38', 'X44', 'X45', 'X46', 'X47', 'X48', 'X54', 'X55', 'X56', 'X57', 'X58', 'X59', 'X60', 'X73')

print("Getting raw clade/isotype freqs")
clade_iso_freqs = identities %>%
  filter(quality & (!restrict | isotype == "iMet")) %>%
  select(match(c('clade', 'isotype', positions), colnames(identities))) %>%
  gather(positions, bases, -clade, -isotype) %>%
  group_by(clade, isotype, positions, bases) %>%
  tally() %>%
  group_by(clade, isotype, positions) %>%
  mutate(freq=n) %>%
  group_by(clade, isotype, positions) %>%
  summarize(A = sum(freq[bases == "A"]),
            C = sum(freq[bases == "C"]),
            G = sum(freq[bases == "G"]),
            U = sum(freq[bases == "U"]),
            Deletion = sum(freq[bases %in% c("-", ".")]), 
            Purine = sum(freq[bases %in% c("A", "G")]),
            Pyrimidine = sum(freq[bases %in% c("C", "U")]),
            Weak = sum(freq[bases %in% c("A", "U")]),
            Strong = sum(freq[bases %in% c("G", "C")]),
            Amino = sum(freq[bases %in% c("A", "C")]),
            Keto = sum(freq[bases %in% c("G", "U")]),
            B = sum(freq[bases %in% c("C", "G", "U")]),
            D = sum(freq[bases %in% c("A", "G", "U")]),
            H = sum(freq[bases %in% c("A", "C", "U")]),
            V = sum(freq[bases %in% c("A", "C", "G")]),
            D = sum(freq[bases %in% c("A", "G", "U")]),
            GC = sum(freq[bases == "G:C"]),
            AU = sum(freq[bases == "A:U"]),
            UA = sum(freq[bases == "U:A"]),
            CG = sum(freq[bases == "C:G"]),
            GU = sum(freq[bases == "G:U"]),
            UG = sum(freq[bases == "U:G"]),
            PairDeletion = sum(freq[bases == "-:-"]), 
            PurinePyrimidine = sum(freq[bases %in% c("A:U", "G:C")]),
            PyrimidinePurine = sum(freq[bases %in% c("U:A", "C:G")]),
            StrongPair = sum(freq[bases %in% c("G:C", "C:G")]),
            WeakPair = sum(freq[bases %in% c("A:U", "U:A")]),
            Wobble = sum(freq[bases %in% c("G:U", "U:G")]),
            Paired = sum(freq[bases %in% c("A:U", "U:A", "C:G", "G:C", "G:U", "U:G")]),
            Bulge = sum(freq[bases %in% c("A:-", "U:-", "C:-", "G:-", "-:A", "-:G", "-:C", "-:U")]),
            Mismatched = sum(freq[bases %in% c("A:A", "G:G", "C:C", "U:U", "A:G", "A:C", "C:A", "C:U", "G:A", "U:C")])
            ) %>%
  mutate(total = A + B + Deletion + Paired + Mismatched + Bulge + PairDeletion) %>%
  melt(id.vars=c("clade", "isotype", "positions", "total")) %>%
  mutate(freq=value/total)

get_isotype_IDE = function(isotype, position, codes) {
  isotype = unique(isotype)
  position = unique(position)
  valid_codes = resolve_code(codes) # returns a vector of all possible combinations of bases. Note that this is not limited to the basic combinations given by the codes variable.
  if (length(isotype) != 1) stop("Multiple isotypes passed to function")
  if (length(position) != 1) stop("Multiple positions passed to function")
  if (length(codes) != 7) return("N/A") # make sure that each clade is represented
  if (length(valid_codes) == 0) return("N/A")
  codes = list(A = "A", C = "C", G = "G", U = "U", Deletion = c("-", "."), Purine = c("A", "G"), Pyrimidine = c("C", "U"), Weak = c("A", "U"), Strong = c("G", "C"), Amino = c("A", "C"), Keto = c("G", "U"), B = c("C", "G", "U"), D = c("A", "G", "U"), H = c("A", "C", "U"), V = c("A", "C", "G"), D = c("A", "G", "U"), GC =  "G:C", AU =  "A:U", UA =  "U:A", CG =  "C:G", GU =  "G:U", UG =  "U:G", PairDeletion =  "-:-", PurinePyrimidine = c("A:U", "G:C"), PyrimidinePurine = c("U:A", "C:G"), StrongPair = c("G:C", "C:G"), WeakPair = c("A:U", "U:A"), Wobble = c("G:U", "U:G"), Paired = c("A:U", "U:A", "C:G", "G:C", "G:U", "U:G"), Bulge = c("A:-", "U:-", "C:-", "G:-", "-:A", "-:G", "-:C", "-:U"), Mismatched = c("A:A", "G:G", "C:C", "U:U", "A:G", "A:C", "C:A", "C:U", "G:A", "U:C")) 
  # for each possible code, check each species
  # this ensures that we consider other IDEs if the most specific one fails the species check. For example, if "Purine" fails, we also consider "V"
  for (code in valid_codes) {
    codes_str = paste0("c('", paste0(codes[[code]], collapse="', '"), "')")
    df = identities %>% 
           select_('species', 'isotype', position) %>%    
           filter_(paste0("isotype == '", isotype, "'")) %>%
           group_by_('species') %>% 
           summarize_(match = paste0("sum(", position, " %in% ", codes_str, ")"),
                      miss = paste0("sum(!(", position, " %in% ", codes_str, "))")) %>%
           mutate(ubiquitous = match > miss)
    if (all(df$ubiquitous)) return(code)
  }
  return("N/A")
}

get_clade_isotype_IDE = function(clade, isotype, position, codes) {
  clade = unique(clade)
  isotype = unique(isotype)
  position = unique(position)
  valid_codes = resolve_code(codes) # returns a vector of all possible combinations of bases. Note that this is not limited to the basic combinations given by the codes variable.
  if (length(clade) != 1) stop("Multiple clades passed to function")
  if (length(isotype) != 1) stop("Multiple clades passed to function")
  if (length(position) != 1) stop("Multiple positions passed to function")
  if (length(valid_codes) == 0) return("N/A")
  codes = list(A = "A", C = "C", G = "G", U = "U", Deletion = c("-", "."), Purine = c("A", "G"), Pyrimidine = c("C", "U"), Weak = c("A", "U"), Strong = c("G", "C"), Amino = c("A", "C"), Keto = c("G", "U"), B = c("C", "G", "U"), D = c("A", "G", "U"), H = c("A", "C", "U"), V = c("A", "C", "G"), D = c("A", "G", "U"), GC =  "G:C", AU =  "A:U", UA =  "U:A", CG =  "C:G", GU =  "G:U", UG =  "U:G", PairDeletion =  "-:-", PurinePyrimidine = c("A:U", "G:C"), PyrimidinePurine = c("U:A", "C:G"), StrongPair = c("G:C", "C:G"), WeakPair = c("A:U", "U:A"), Wobble = c("G:U", "U:G"), Paired = c("A:U", "U:A", "C:G", "G:C", "G:U", "U:G"), Bulge = c("A:-", "U:-", "C:-", "G:-", "-:A", "-:G", "-:C", "-:U"), Mismatched = c("A:A", "G:G", "C:C", "U:U", "A:G", "A:C", "C:A", "C:U", "G:A", "U:C")) 
  # for each possible code, check each species
  # this ensures that we consider other IDEs if the most specific one fails the species check. For example, if "Purine" fails, we also consider "V"
  for (code in valid_codes) {
    codes_str = paste0("c('", paste0(codes[[code]], collapse="', '"), "')")
    df = identities %>% 
           select_('species', 'isotype', 'clade', position) %>%    
           filter_(paste0("clade == '", clade, "' & isotype == '", isotype, "'")) %>%
           group_by_('species') %>% 
           summarize_(match = paste0("sum(", position, " %in% ", codes_str, ")"),
                      miss = paste0("sum(!(", position, " %in% ", codes_str, "))")) %>%
           mutate(ubiquitous = match > miss)
    if (all(df$ubiquitous)) return(code)
  }
  return("N/A")
}

df = data.frame()
for (cutoff in c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99)) {
  print(paste0("Processing freqs at a cutoff of ", cutoff))
  best_freqs = clade_iso_freqs %>%
    group_by(isotype, positions, clade, variable) %>% # remove duplicates
    summarize(count=sum(value), freq=sum(value)/sum(total)) %>%
    filter(freq > cutoff) %>%
    mutate(cutoff=cutoff) %>%
    group_by(isotype, clade, positions) %>%
    filter(row_number(freq) == 1)
  isotype_specific = best_freqs %>%
    group_by(isotype, positions, cutoff) %>% 
    summarize(identity = get_isotype_IDE(isotype, positions, variable)) %>%
    mutate(clade='Eukarya') %>%
    filter(identity != "N/A")
  clade_isotype_specific = best_freqs %>%
    group_by(clade, isotype, positions, cutoff) %>% 
    summarize(identity = get_clade_isotype_IDE(clade, isotype, positions, variable)) %>%
    filter(identity != 'N/A')

  if (nrow(df) == 0) df = rbind(isotype_specific, clade_isotype_specific)
  else df = rbind(df, isotype_specific, clade_isotype_specific)
}
write.table(df, file='freqs.tsv', sep='\t', quote=FALSE, row.names=FALSE)
