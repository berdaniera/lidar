publicHeaderDescription <- function() {
    hd <- structure(list(Item = c("File Signature (\"LASF\")",
"(1.1) File Source ID", "(1.1) Global Encoding",
"(1.1) Project ID - GUID data 1", "(1.1) Project ID - GUID data 2",
"(1.1) Project ID - GUID data 3", "(1.1) Project ID - GUID data 4",
"Version Major", "Version Minor", "(1.1) System Identifier",
"Generating Software", "(1.1) File Creation Day of Year",
"(1.1) File Creation Year", "Header Size", "Offset to point data",
"Number of variable length records",
"Point Data Format ID (0-99 for spec)", "Point Data Record Length",
"Number of point records", "Number of points by return",
"X scale factor", "Y scale factor", "Z scale factor", "X offset",
"Y offset", "Z offset", "Max X", "Min X", "Max Y", "Min Y", "Max Z",
"Min Z"), Format = c("char[4]", "unsigned short", "unsigned short",
"unsigned long", "unsigned short", "unsigned short",
"unsigned char[8]", "unsigned char", "unsigned char", "char[32]",
"char[32]", "unsigned short", "unsigned short", "unsigned short",
"unsigned long", "unsigned long", "unsigned char", "unsigned short",
"unsigned long", "unsigned long[5]", "double", "double", "double",
"double", "double", "double", "double", "double", "double", "double",
"double", "double"), Size = c("4 bytes", "2 bytes", "2 bytes",
"4 bytes", "2 byte", "2 byte", "8 bytes", "1 byte", "1 byte",
"32 bytes", "32 bytes", "2 bytes", "2 bytes", "2 bytes", "4 bytes",
"4 bytes", "1 byte", "2 bytes", "4 bytes", "20 bytes", "8 bytes",
"8 bytes", "8 bytes", "8 bytes", "8 bytes", "8 bytes", "8 bytes",
"8 bytes", "8 bytes", "8 bytes", "8 bytes", "8 bytes"), Required =
c("*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*",
"*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*", "*",
"*", "*", "*", "*", "*")), .Names = c("Item", "Format", "Size",
"Required"), row.names = 2:33, class = "data.frame")
    hd$what <- ""
    hd$what[grep("unsigned", hd$Format)] <- "integer"
    hd$what[grep("char", hd$Format)] <- "raw"
    hd$what[grep("short", hd$Format)] <- "integer"
    hd$what[grep("long", hd$Format)] <- "integer"
    hd$what[grep("double", hd$Format)] <- "numeric"
    hd$signed <- TRUE
    hd$signed[grep("unsigned", hd$Format)] <- FALSE
    ## number of values in record
    hd$n <- as.numeric(gsub("[[:alpha:][:punct:]]", "", hd$Format))
    hd$n[hd$what == "character"] <- 1
    hd$n[is.na(hd$n)] <- 1
    ## size of record
    hd$Hsize <- as.numeric(gsub("[[:alpha:]]", "", hd$Size))
    ## size of each value in record
    hd$Rsize <- hd$Hsize / hd$n
    hd$Rsize[hd$what == "raw"] <- 1
    hd$n[hd$what == "raw"] <- hd$Hsize[hd$what == "raw"]
    hd
}

readLAS <-
function(lasfile, skip = 0, nrows = NULL, returnSP = FALSE, returnHeaderOnly = FALSE) {
    ## TODO:
    ## generalize to any version
    ## figure out what this gpstime is to convert to POSIXct .. .
    ## how do we get coordinate system?
    ## bits after Intensity
    ## can we be more efficient to "seek" without actually reading bytes on windows?

    ## DONE
    ## ensure the entire header is read, from Header Size
    ## parse header
    ## provide chunked read

    hd <- publicHeaderDescription()
    pheader <- vector("list", nrow(hd))
    names(pheader) <- hd$Item
    con <- file(lasfile, open = "rb")
    isLASFbytes <- readBin(con, "raw", size = 1, n = 4, endian = "little")
    pheader[[hd$Item[1]]] <- readBin(isLASFbytes, "character", size = 4, endian = "little")
    if (! pheader[[hd$Item[1]]] == "LASF") {
        stop("Not a valid LAS file")
    }
    for (i in 2:nrow(hd)) {
        pheader[[hd$Item[i]]] <- readBin(con, what = hd$what[i], signed = hd$signed[i], size = hd$Rsize[i], endian = "little", n = hd$n[i])
        #print(names(pheader)[i])
        #print(pheader[[hd$Item[i]]])
    }
    close(con)
    ## read the data
    numberPointRecords <- pheader[["Number of point records"]]
    offsetToPointData <- pheader[["Offset to point data"]]
    pointDataRecordLength <-pheader[["Point Data Record Length"]]
    xyzScaleOffset <- cbind(unlist(pheader[c("X scale factor", "Y scale factor", "Z scale factor")]),
                            unlist(pheader[c("X offset", "Y offset", "Z offset")]))


    if (returnHeaderOnly) return(pheader)

    ## totaln <- pheader$`Number of point records`
    ## mm <- matrix(raw(), nrow = totaln, ncol = pointDataRecordLength)


    con <- file(lasfile, open = "rb")
    junk <- readBin(con, "raw", size = 1, n = offsetToPointData)

    ## deal with rows to skip and max rows to be read
    if (skip > 0) {
        ## seek is unreliable on windows, or I'm using it incorrectly
        ## so we junk the bytes to skip
        junk <- readBin(con, "raw", size = 1, n = pointDataRecordLength * skip)
	numberPointRecords <- numberPointRecords - skip
        #pos <- seek(con, where = pointDataRecordLength * skip)
        # print(c(pos = seek(con), skip = skip, where = pointDataRecordLength * skip))
    }
    if (!is.null(nrows)) {
	if (numberPointRecords > nrows) numberPointRecords <- nrows
    }

    if (numberPointRecords < 1) stop("no records left to read")

    allbytes <- matrix(readBin(con, "raw", n = pointDataRecordLength * numberPointRecords, size = 1, endian = "little"),
                           ncol= pointDataRecordLength, nrow = numberPointRecords, byrow = TRUE)
    close(con)
    mm <- matrix(readBin(t(allbytes[,1:(3*4)]), "integer", size = 4, n = 3 * numberPointRecords, endian = "little"), ncol = 3, byrow = TRUE)
     gpstime <- NULL
    if (ncol(allbytes) == 28) gpstime <- readBin(t(allbytes[ , 21:28]), "numeric", size = 8, n = numberPointRecords, endian = "little")

    intensity <- readBin(t(allbytes[, 13:14]), "integer", size = 2, n = numberPointRecords, signed = FALSE, endian = "little")
    mm[,1] <- mm[ ,1] * xyzScaleOffset[1,1] + xyzScaleOffset[1, 2]
    mm[,2] <- mm[ ,2] * xyzScaleOffset[2,1] + xyzScaleOffset[2, 2]
    mm[,3] <- mm[ ,3] * xyzScaleOffset[3,1] + xyzScaleOffset[3, 2]
    colnames(mm) <- c("x", "y", "z")
    if (returnSP) {
        require(sp)
        SpatialPoints(cbind(mm, gpstime, intensity))
    } else {
        cbind(mm, gpstime, intensity)
    }
}
