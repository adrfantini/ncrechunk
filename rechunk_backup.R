#!/usr/bin/env Rscript

#May 2018
#afantini@ictp.it

suppressPackageStartupMessages(library(ncdf4))
suppressPackageStartupMessages(library(optparse))

tellme = function(...) {
    message(" ::: ", format(Sys.time(), "%H:%M:%S"), " ::: ", ...)
}

option_list = list(make_option(c("-m", "--memfree"),
                                type="integer",
                                default=1024,
                                help="Amount of memory to keep free, in MB [default: %default]"),
                    make_option(c("-s", "--spacechunk"),
                                type="integer",
                                default=20,
                                help="The size of the non-time chunking, the same for all non-time dimensions. A higher value will have much faster conversion speed and better compression, but will be slower when reading a single timeseries [default: %default]"),
                    make_option(c("-v", "--variable"),
                                type="character",
                                help="Variable to rechunk"),
                    make_option("--timevar",
                                default="time",
                                type="character",
                                help="Name of the time variable [default %default]"),
                    make_option("--suffix",
                                type="character",
                                default="_rechunk",
                                help="String to attach as suffix to the variable name [default %default]"),
                    make_option("--maxchunk",
                                type="integer",
                                default=365*24,
                                help="Maximum size of time chunk [default %default]")
)
parser = OptionParser(
    usage = "%prog -v VARIABLE FILENAME",
    option_list=option_list,
    epilogue=paste0("This script rechunks a given netcdf variable from slice to timeseries. Works around memory limitations by splitting into several jobs. ASSUMES SIGNED VALUES! For feature requests and bug reports, please write to:
    afantini@ictp.it")
)

## TODO Set a maxmem argument, instead of memfree

#Gather input arguments
arguments = parse_args(parser, positional_arguments = 1)
opt = arguments$options
fn = arguments$args
membuffer = opt$memfree
suffix = opt$suffix
maxchunk = opt$maxchunk
var = opt$variable
spacechunk = opt$spacechunk
timevar = opt$timevar

if (is.null(var)) {print_help(parser) ; stop("-v argument is mandatory")}

tellme("Opening input file ", fn)
nc = nc_open(fn, write=TRUE)
if (! var %in% names(nc$var)) stop("The specified variable does not exist in the input file")

#get the wanted chunk specifications
chunknames = sapply(nc$var[[var]]$dim, `[[`, "name")

bytes = switch(nc$var[[var]]$prec,
    "short" = 4,
    "byte" = 4,
    "integer" = 4,
    "float" = 8,
    "double" = 8,
    8
)
overhead = 6
# memfree = round(as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern=TRUE)) / 1024) +
memfree = round(as.numeric(system("awk '/MemAvailable/ {print $2}' /proc/meminfo", intern=TRUE)) / 1024)
memreqd = round(prod(nc$var[[var]]$size) * bytes * overhead / 1024**2)
#overhead is an empirical constant to estimate the amount of stupidity that is going on in this code / in netcdf's library

if (memfree < memreqd + membuffer) {
    nj = memreqd %/% (memfree - membuffer) + 1
    if (nj < 1) stop("ERROR! Not enough memory to do anything!")
    tellme("Not enough free RAM to perform the rechunking in one single job: ", memfree, "MB free, ", memreqd + membuffer, "MB requested (", membuffer, "MB kept as free buffer)
    Will create ", nj, " separate time chunks")
} else {
    tellme("Starting. Expected RAM usage: ~ ", memreqd, " MB")
    nj = 1
}

cs = rep(spacechunk, length(chunknames))
nT = nc$dim[[timevar]]$len
whichT = which(chunknames==timevar)
Tchunk = ceiling(nT / nj)
cs[whichT] = Tchunk

start = seq(1, nT, by = Tchunk)
count = rep(Tchunk, length(start))
count[length(count)] = nT - start[length(count)] + 1

tellme("Time chunks used:")
tsdf = data.frame(start=start, end = start+count-1, length=count)
if (dim(tsdf)[1] > 20) {
    print(head(tsdf, 10))
    cat("...\n")
    print(tail(tsdf, 10))
} else {
    print(tsdf)
}

tellme("Waiting 2 minutes to start, to give you time to reflect on your sins")
Sys.sleep(120)

newname = paste0(var, suffix)
if (newname %in% names(nc$var)) {
    renamed = paste0(newname, "_old")
    while (renamed %in% names(nc$var)) renamed = paste0(renamed, "_old") #veeeeeeery old?
    tellme("The file already has a variable named ", newname, ", renaming it to ", renamed)
    nc = ncvar_rename(nc, newname, renamed)
}

newvar = ncvar_def(
    name = newname,
    units = nc$var[[var]]$units,
    dim = nc$var[[var]]$dim,
    prec= nc$var[[var]]$prec,
    missval = nc$var[[var]]$missval,
    longname = nc$var[[var]]$longname,
    compression = 1,
    chunksizes = cs
    )
nc = ncvar_add(nc, newvar)

invisible(gc())

for (i in 1:nj) {
    istart = rep(1, length(cs))
    icount = rep(-1, length(cs))
    istart[whichT] = start[i]
    icount[whichT] = count[i]
    tellme("Writing chunk ", i, " of ", nj)
    ncvar_put(nc, newvar,
              start = istart,
              count = icount,
              vals = ncvar_get(nc, var, collapse_degen=FALSE, raw_datavals=TRUE,
                               start = istart,
                               count = icount
                               )
              )

}

nc_close(nc)

invisible(gc())

nc = nc_open(fn, write=TRUE)
atts <- ncatt_get(nc, var)
if (length(atts) > 0) {
    invisible(lapply(names(atts), function(n) {
            try( ncatt_put( nc, newvar, n, atts[[n]] ) )
    }))
}
## add global var, history
nc_close(nc)

invisible(gc())

tellme("All done!")
