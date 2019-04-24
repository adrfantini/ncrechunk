#!/usr/bin/env Rscript

#May 2018
#afantini@ictp.it

suppressPackageStartupMessages(library(ncdf4))
suppressPackageStartupMessages(library(optparse))

option_list = list(make_option(c("-m", "--memfree"),
                                type="integer",
                                default=1024,
                                help="Amount of memory to keep free, in MB [default: %default]"),
                    make_option(c("-c", "--chunkspec"),
                                type="character",
                                help="Comma separated list of chunk sizes. Values <=0 mean maximum chunking along that dimension (e.g. lat/-1,lon/-1,time/1)"),
                    make_option(c("-v", "--variable"),
                                type="character",
                                help="Variable to rechunk"),
                    make_option("--suffix",
                                type="character",
                                default="_rechunk",
                                help="String to attach as suffix to the variable name [default %default]")
)
parser = OptionParser(
    usage = "%prog -v VARIABLE -c CHUNKSPEC FILENAME",
    option_list=option_list,
    epilogue=paste0("This script rechunks a given netcdf variable in a file according to the given chunk specifications. ASSUMES SIGNED VALUES! For feature requests and bug reports, please write to:
    afantini@ictp.it")
)
#Gather input arguments
arguments = parse_args(parser, positional_arguments = 1)
opt = arguments$options
fn = arguments$args
membuffer = opt$memfree
chunkspec = opt$chunkspec
suffix = opt$suffix
var = opt$variable

if (is.null(var)) {print_help(parser) ; stop("-v argument is mandatory")}
if (is.null(chunkspec)) {print_help(parser) ; stop("-c argument is mandatory")}

message("Opening input file ", fn)
nc = nc_open(fn, write=TRUE)
if (! var %in% names(nc$var)) stop("The specified variable does not exist in the input file")

#get the wanted chunk specifications
cs = strsplit(chunkspec, ",")[[1]]
cs = as.data.frame(do.call("cbind",
                            lapply(cs, function(str) unlist(strsplit(str, "/")) )
                    ), stringsAsFactors = FALSE)
colnames(cs) = cs[1,]
cs = cs[-1,]
cs[] = as.integer(cs)

for (n in which(cs <= 0)) { #check if values are <=0, and substitute for full chunk
    varn = colnames(cs)[n]
    cs[,n] = nc$dim[[varn]]$len
}

chunknames = sapply(nc$var[[var]]$dim, `[[`, "name")

bytes = switch(nc$var[[var]]$prec,
    "short" = 4,
    "byte" = 4,
    "integer" = 4,
    "float" = 8,
    "double" = 8,
    8
)
overhead = 7
memfree = round(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE)) / 1024)
memreqd = round(prod(nc$var[[var]]$size) * bytes * overhead / 1024**2)
#overhead is an empirical constant to estimate the amount of stupidity that is going on in this code
    
if (memfree < memreqd + membuffer) {
    stop("Not enough free RAM to perform the rechunking: ", memfree, "MB free, ", round(memreqd + membuffer), "MB requested")
} else {
    message("Starting. Expected RAM usage: ~ ", memreqd, " MB")
}

newvar = ncvar_def(
    name = paste0(var, suffix) ,
    units = nc$var[[var]]$units,
    dim = nc$var[[var]]$dim,
    prec= nc$var[[var]]$prec,
    missval = nc$var[[var]]$missval,
    longname = nc$var[[var]]$longname,
    compression = 1,
    chunksizes = as.integer(cs[,chunknames])
    )
nc = ncvar_add(nc, newvar)

tmp = ncvar_get(nc, var, collapse_degen=FALSE, raw_datavals=TRUE)
message("Read everything. Writing...")
ncvar_put(nc, newvar, vals = tmp)

nc_close(nc)
nc = nc_open(fn, write=TRUE)
atts <- ncatt_get(nc, var)
if (length(atts) > 0) {
    invisible(lapply(names(atts), function(i) {
            try(ncatt_put(nc, newvar, i, atts[[i]]))
    }))
}

nc_close(nc)

message("All done!")
