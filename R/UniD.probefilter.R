#' Probe and sample filter
#'
#' @description Used to filter out probes with suspect quality
#' issues for following categories: (1) probes located on ChrX/Y; (2) probes
#' may affected by SNP; (3) probes may mapped to multiple locations; (4)
#' probes not targeted CpG sites; (5) probes on 450k platform but not
#' available on EPIC platform; (6) probes with large proportion of missing
#' values.
#' It can also filter out samples with large proportion of missing
#' values. Missing values may caused by non-significant detection p-value
#' or less beadcount per probe.
#'
#' @param Beta.raw data frame generated from the \code{UniD.dataqc()}
#' @param filterXY whether filter out probes located on Chromosome
#' X and Y. Default value is TRUE.
#' @param filterSNPHit whether filter out probes may affected by SNP.
#' Default is TRUE. Probe list adapted from Zhou W, Nucleic Acids Research,
#'  2017.
#' @param filterMultiHit whether filter out probes can mapped to
#' multi locations. Default is TRUE. Probe list adapted from Nordlund J,
#'  Genome Biology, 2013.
#' @param filterNonCG whether filter out probes which are not targeting
#' CpG sites. Default is FALSE.
#' @param filterNonEpic whether filter outprobes which are available on
#' 450k platform but not available on EPIC platform. Default is T. Highly
#' recommended if building models with data from 450k platform.
#' @param filterSample whether filter out samples with high proportion
#' of missing values. Default is TRUE.
#' @param filterSample.cut if \code{filterSample} = T, the threshold for
#' high proportion of missing values per sample. Default is 0.1.
#' @param filterProbe whether filter out probes with high proportion
#' of missing values. Default is TRUE. Carefully usage with small number of
#' sample size.
#' @param filterProbe.cut If \code{filterProbe} = T, the threshold for
#' high proporition of missing values per probe. Default is 0.05.
#' @inheritParams UniD.load
#' @return A data frame with probes after filtering in Beta value
#' filtering
#' @examples
#' \dontrun{
#' Beta.clean <- UniD.probefilter(Beta.raw, outDir = "~/Desktop/output/",
#' filterXY = T, filterSNPHit = T, filterMultiHit = T, filterNonCG = F,
#' filterNonEpic = T, arrayType = "450k", filterSample = T, filterSample.cut
#' = 0.1, filterProbe = F, write = T)
#'
#' Beta.clean <- UniD.probefilter(Beta.raw, outDir = NULL,
#' filterXY = T, filterSNPHit = T, filterMultiHit = T, filterNonCG = F,
#' filterNonEpic = F, arrayType = "EPIC", filterSample = T, filterSample.cut
#' = 0.1, filterProbe = F, write = F)
#' }
#' @references
#' Zhou, W., et al. (2017). "Comprehensive characterization, annotation and innovative use of Infinium DNA #' methylation BeadChip probes." Nucleic Acids Res 45(4): e22.
#'
#' Nordlund, J., et al. (2013). "Genome-wide signatures of differential DNA methylation in pediatric acute #' lymphoblastic leukemia." Genome Biol 14(9): r105.
#' @import wateRmelon

#' @export
UniD.probefilter <- function(Beta.raw,
                             outDir,
                             filterXY = TRUE,
                             filterSNPHit = TRUE,
                             filterMultiHit = TRUE,
                             filterNonCG = FALSE,
                             filterNonEpic = TRUE,
                             arrayType = c("450k", "EPIC"),
                             filterSample = TRUE,
                             filterSample.cut = 0.1,
                             filterProbe = FALSE,
                             filterProbe.cut = 0.05,
                             write)
{
  message("")
  message("====Probe Filter Start====")

  if(exists("Beta.raw") == FALSE){
    stop("Beta.raw object is missing")
  }

  #delete chrXY
  if(filterXY){
    if(arrayType == "EPIC"){
      Beta.raw2 <- Beta.raw[which(! rownames(Beta.raw) %in% Epic_chrXY),]
    }else if(arrayType == "450k"){
      Beta.raw2 <- Beta.raw[which(! rownames(Beta.raw) %in% k450_chrXY),]
    }
    message(paste0("Delete ", nrow(Beta.raw) - nrow(Beta.raw2),
                   " probes on ChrX/Y for ", arrayType, " Platform."))
    Beta.raw <- Beta.raw2
    rm(Beta.raw2)
  }

  #delete SNP mask
  if(filterSNPHit){
    if(arrayType == "EPIC"){
      Beta.raw2 <- Beta.raw[which(! rownames(Beta.raw) %in% Epic_snp),]
    }else if(arrayType == "450k"){
      Beta.raw2 <- Beta.raw[which(! rownames(Beta.raw) %in% k450_snp),]
    }
    message(paste0("Delete ", nrow(Beta.raw) - nrow(Beta.raw2),
                   " probes affected by SNP for ",
                   arrayType, " Platform. (Zhou W, Nucleic Acids Research, 2017)"))
    Beta.raw <- Beta.raw2
    rm(Beta.raw2)
  }

  #delete multi-hit(450k version)
  if(filterMultiHit){
    Beta.raw2 <- Beta.raw[which(! rownames(Beta.raw) %in% multi.hit),]
    message(paste0("Delete ", nrow(Beta.raw) - nrow(Beta.raw2),
                   " probes mapped to multiple site for ",
                   arrayType, " Platform. (Nordlund J, Genome Biology, 2013)"))
    Beta.raw <- Beta.raw2
    rm(Beta.raw2)
  }

  #delete probes in 450k but not in EPIC
  if(filterNonEpic){
    if(arrayType == "450k"){
      Beta.raw2 <- Beta.raw[which(! rownames(Beta.raw) %in% NonEpic),]
      message(paste0("Delete ", nrow(Beta.raw) - nrow(Beta.raw2),
                     " probes not available on Epic platform for ",
                     arrayType, " Platform"))
      Beta.raw <- Beta.raw2
      rm(Beta.raw2)
    }
  }

  #delete probes not targeted for CG sites
  if(filterNonCG){
    Beta.raw2 <- Beta.raw[which(! grepl("ch.", row.names(Beta.raw))),]
    message(paste0("Delete ", nrow(Beta.raw) - nrow(Beta.raw2),
                   " probes target non-CG site for ", arrayType, " Platform"))
    Beta.raw <- Beta.raw2
    rm(Beta.raw2)
  }

  #filter samples with too many NAs among probes
  if(filterSample){
    p.fail <- data.frame(Fail.Frac.NA = colSums(is.na(Beta.raw))/nrow(Beta.raw))
    print("The failed fraction per sample (after all probe filters): ")
    print(p.fail)

    filter.ID <- which(p.fail$Fail.Frac.NA > filterSample.cut)

    if(length(filter.ID) > 0){
      Beta.raw2 <- Beta.raw[, - filter.ID]
      message(paste0("Delete sample ", paste(colnames(Beta.raw)[filter.ID],
                                             collapse = ";"),
                     " due to missing value > ", filterSample.cut,
                     " among probes"))
      Beta.raw <- Beta.raw2
      rm(Beta.raw2)
    }
  }


  #filter probes with too many NAs among samples
  if(filterProbe){
    s.fail <- data.frame(Fail.Frac.NA.perSample = rowSums(is.na(Beta.raw))
                         /ncol(Beta.raw))
    print("The failed fraction per probe (after all probe filters): ")
    print(summary(s.fail))

    filter.ID <- which(s.fail$Fail.Frac.NA.perSample > filterProbe.cut)
    if(length(filter.ID) > 0){
      Beta.raw2 <- Beta.raw[ - filter.ID , ]
      message(paste0("Delete ", nrow(Beta.raw) - nrow(Beta.raw2),
                     " probe due to missing value > ", filterProbe.cut,
                     " among samples"))
      Beta.raw <- Beta.raw2
      rm(Beta.raw2)
    }
  }

  message("")
  message("Final data matrix is: ", ncol(Beta.raw), " samples * ",
          nrow(Beta.raw), " probes.")

  if(write == TRUE){
    Beta.clean <- Beta.raw
    save(Beta.clean, file = paste0(outDir,"/UniD_Beta.clean.RData"))
  }

  message("====Probe Filter Finish====")
  message("")
  return(Beta.clean)
}
