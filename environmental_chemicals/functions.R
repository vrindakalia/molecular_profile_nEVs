replaceLODzero <- function(x) "[<-"(x, x == "<LOD", "0")

replaceNAzero <- function(x) "[<-"(x, is.na(x), "0")

replaceNFNA <- function(x) "[<-"(x, x=="N/F", "NA")

replacezero <- function(x) "[<-"(x, !x|x==0, min(x[x>0], na.rm = T)/2)

codeLODNF <- function(x) {
    y = case_when(x == "<LOD" ~ 0,
                  x == "N/F" ~ 0,
                  x != "<LOD|N/F" ~1,
                  is.na(x) ~ 0)
    return(y)
}

FindMin <- function(x) (min(x[x>0]))
