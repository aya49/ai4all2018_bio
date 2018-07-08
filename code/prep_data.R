# title: "perpare sequence and metadata"
# author: "aya43@sfu.ca"
# date: '2018-07-07'

# load packages
pkgs = c("ape", "data.table")
pkgs_ui = setdiff(pkgs, rownames(installed.packages()))
if (length(pkgs_ui) > 0) install.packages(pkgs_ui, verbose=F)
sapply(pkgs, require, character.only=TRUE)

# set paths (change "root" to your directory)
root = "~/projects/ai4all-sfu_bio"

seq_dir_ = paste(root, "/data_/FASTA.fa", sep="") # original genetic sequence
meta_dir_ = paste(root, "/data_/meta.tsv", sep="") # original metadata

seq_dir = paste(root, "/data/FASTA.fa", sep="")
meta_dir = paste(root, "/data/meta.csv", sep="")

result_dir = paste(root, "/result", sep=""); dir.create(result_dir, showWarnings=F)
align_dir = paste(result_dir, "/alignment_mafft.fa", sep="")

# load data
sequ_ = read.dna(seq_dir_, format="fasta")
meta_ = fread(meta_dir_, data.table=F)

strains_ = sapply(strsplit(names(sequ_),"[-]"), function(x) x[1])
dates_ = as.character(as.Date(sapply(strsplit(names(sequ_),"[-]"), function(x) x[2])))

# order data
ord = order(dates_, decreasing=TRUE)
sequ = sequ_[ord]
strains = strains_[ord]
dates = dates_[ord]
dates_split = strsplit(dates,"[-]")
years = sapply(dates_split, function(x) x[1])
months = sapply(dates_split, function(x) x[2])
days = sapply(dates_split, function(x) x[3])

# clean metadata
meta = meta_[match(strains, meta_[,"Sequence Accession"]),]
meta[grepl("-N/A-", meta[,"Flu Season"]), "Flu Season"] = NA
current_year = as.numeric(substr(as.character(Sys.Date()), 3, 4))
meta[!is.na(meta[,"Flu Season"]), "Flu Season"] = 
  sapply(strsplit(meta[!is.na(meta[,"Flu Season"]),"Flu Season"],"[-]"), function(x) {
    xn = as.numeric(x)
    return( paste(ifelse(xn[1]<=current_year, "20", "19"), x[1], "-", 
                  ifelse(xn[2]<=current_year, "20", "19"), x[2], sep="") )
  })
meta[is.na(meta[,"Flu Season"]), "Flu Season"] = 
  mapply(function(y,x) {
    y = as.numeric(y)
    yp = ifelse(as.numeric(x)>6, 0, -1)
    return( paste(y+yp, "-", y+yp+1, sep="") )
  }, years[is.na(meta[,"Flu Season"])], months[is.na(meta[,"Flu Season"])])

meta = meta[,c("Sequence Accession", "Subtype", "Collection Date", "Flu Season", "Country", "State/Province")] #only keep a few columns
colnames(meta) = c("strain", "subtype", "date", "season", "country", "state")

# save processed data
write.FASTA(sequ, file=seq_dir)
write.csv(meta, file=meta_dir)





