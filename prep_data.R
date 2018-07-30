# title: "perpare sequence and metadata"
# author: "aya43@sfu.ca"
# date: '2018-07-07'

# load packages
pkgs = c("ape", "data.table", "RJSONIO")
pkgs_ui = setdiff(pkgs, rownames(installed.packages()))
if (length(pkgs_ui) > 0) install.packages(pkgs_ui, verbose=F)
sapply(pkgs, require, character.only=TRUE)
options(geonamesUsername="silverloc123")

# set paths (change "root" to your directory)
root = "~/projects/ai4all-sfu_bio"

seq_dir_ = paste(root, "/data_/FASTA.fa", sep="") # original genetic sequence
meta_dir_ = paste(root, "/data_/meta_fludb.tsv", sep="") # original metadata
meta2_dir_ = paste(root, "/data_/meta_ncbi.txt", sep="") # original metadata

seq_dir = paste(root, "/data/FASTA.fa", sep="")
meta_dir = paste(root, "/data/meta.csv", sep="")
align_dir = paste(root, "data/alignment_mafft.fa", sep="")

# load data
sequ_ = read.dna(seq_dir_, format="fasta")
meta_ = fread(meta_dir_, data.table=F, stringsAsFactors=F)
meta2_ = read.csv(meta2_dir_, stringsAsFactors=F)

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

meta = meta_[,c("Sequence Accession", "Subtype", "Collection Date", "Country", "State/Province", "Flu Season")] #only keep a few columns
colnames(meta) = c("strain", "subtype", "date", "country", "state", "season")
meta = meta[match(strains, meta[,"strain"]),]
meta_m = as.matrix(meta); meta_m[meta_m=="-N/A-"] = NA
meta_m[,"date"] = dates

# patch up metadata
if (sum(is.na(meta_m[,"strain"])) > 0) {
  meta_patch_ = meta2_[match(strains[is.na(meta_m[,"strain"])], meta2_[,"accession"]),]
  meta_patch = cbind(meta_patch_[,c("accession", "serotype", "date", "country")], matrix(NA, nrow=nrow(meta_patch_), ncol=2))
  meta_m[is.na(meta_m[,"strain"]),] = as.matrix(meta_patch)
}

# clean metadata
meta = as.data.frame(meta_m, stringsAsFactors=F)
current_year = as.numeric(substr(as.character(Sys.Date()), 3, 4))
meta[!is.na(meta[,"season"]), "season"] = 
  sapply(dates_split[!is.na(meta[,"season"])], function(x) {
    xn = as.numeric(x)
    return( paste(ifelse(xn[1]<=current_year, "20", "19"), x[1], "-", 
                  ifelse(xn[2]<=current_year, "20", "19"), x[2], sep="") )
  })
meta[is.na(meta[,"season"]), "season"] = 
  mapply(function(y,x) {
    y = as.numeric(y)
    yp = ifelse(as.numeric(x)>6, 0, -1)
    return( paste(y+yp, "-", y+yp+1, sep="") )
  }, years[is.na(meta[,"season"])], months[is.na(meta[,"season"])])

# retrieve geographic coordinates
geo_query = as.character(meta[,"country"])
geo_query[!is.na(meta[,"state"])] = paste(meta[!is.na(meta[,"state"]),"state"],", ", geo_query[!is.na(meta[,"state"])], sep="")

# https://stackoverflow.com/questions/32504880/street-address-to-geolocation-lat-long
address_coord = function(address) {
  url = URLencode(paste("http://maps.google.com/maps/api/geocode/json?address=", address, "&sensor=false", sep=""))
  x = RJSONIO::fromJSON(url, simplify=FALSE)
  Sys.sleep(0.21)  # API only allows 5 requests per second
  while (x$status!="OK") {
    x = RJSONIO::fromJSON(url, simplify=FALSE)
    Sys.sleep(0.21)
  }
  return(unlist(x$results[[1]]$geometry$location))
}

coords = sapply(unique(geo_query), address_coord)
coords_table = t(coords)

meta = cbind(meta, coords_table[geo_query,])
meta[,"date"] = sapply(strsplit(as.character(meta[,"date"]), "[-/]"), function(x) paste(x, collapse="-"))



# save processed data
write.FASTA(sequ, file=seq_dir)
write.csv(meta, file=meta_dir)





