import os
import pandas as pd
import Bio.SeqIO as Seq
import numpy as np
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

def _get_r_parallel_cores():
    # LncFinder parallel backend may hang/fail in some environments (WSL, restricted sockets).
    # Default to single core; allow override via env var.
    value = os.getenv("NCRESIS_R_CORES", "1")
    try:
        cores = int(value)
    except ValueError:
        cores = 1
    return max(1, cores)

def extract_SSfeatures(fastapath):
    parallel_cores = _get_r_parallel_cores()
    robjects.r('''
                extract_SSfeatures <- function(fastapath, parallel_cores){
                  library(LncFinder)
                  demo_DNA.seq <- seqinr::read.fasta(fastapath)
                  Seqs <- LncFinder::run_RNAfold(demo_DNA.seq, RNAfold.path = "RNAfold", parallel.cores = parallel_cores)
                  result_2 <- LncFinder::extract_features(Seqs, label = NULL, SS.features = TRUE,format = "SS", frequencies.file = "human", parallel.cores = parallel_cores)
                  res2 <- result_2[,c(12:19)]
                  return(res2)
                }''')
    sstruc = robjects.r['extract_SSfeatures'](fastapath, parallel_cores)
    sstruc.columns = ["SLDLD: Structural logarithm distance to lncRNA of acguD", "SLDPD: Structural logarithm distance to pcRNA of acguD", "SLDRD: Structural logarithm distance acguD ratio", "SLDLN: Structural logarithm distance to lncRNA of acguACGU", "SLDPN: Structural logarithm distance to pcRNA of acguACGU", "SLDRN: Structural logarithm distance acguACGU ratio","SDMFE: Secondary structural minimum free energy", "SFPUS: Secondary structural UP frequency paired-unpaired"]
    return sstruc


def makeEIIP(fastapath):
    robjects.r('''
                makeEIIP <- function(fastapath){
                  library(LncFinder)
                  demo_DNA.seq <- seqinr::read.fasta(fastapath)
                  result_1 <- compute_EIIP(
                                demo_DNA.seq,
                                label = NULL,
                                spectrum.percent = 0,
                                quantile.probs = seq(0, 1, 0.25)
                                )
                  return(result_1)
                }''')
    sstruc = robjects.r['makeEIIP'](fastapath)
    sstruc.columns = ['EipSP: Electron-ion interaction pseudopotential signal peak','EipAP: Electron-ion interaction pseudopotential average power','EiSNR: Electron-ion interaction pseudopotential signal/noise ratio','EiPS0: Electron-ion interaction pseudopotential spectrum 0','EiPS1: Electron-ion interaction pseudopotential spectrum 0.25','EiPS2: Electron-ion interaction pseudopotential spectrum 0.5','EiPS3: Electron-ion interaction pseudopotential spectrum 0.75','EiPS4: Electron-ion interaction pseudopotential spectrum 1']
    return sstruc

def makeORFEucDist(fastapath):
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri
    current_script_dir = os.path.dirname(os.path.abspath(__file__))
    abs_cds_path = os.path.join(current_script_dir, 'Data', 'gencode.v34.pc_transcripts_test.fa')
    abs_lncrna_path = os.path.join(current_script_dir, 'Data', 'gencode.v34.lncRNA_transcripts_test.fa')

    robjects.r('''
                makeORFEucDist <- function(fastapath, cds_path, lncrna_path){
                  library(LncFinder)

                  cds.seq = seqinr::read.fasta(cds_path)
                  lncRNA.seq = seqinr::read.fasta(lncrna_path)   

                  referFreq <- make_referFreq(
                            cds.seq,
                            lncRNA.seq,
                            k = 6,
                            step = 1,
                            alphabet = c("a", "c", "g", "t"),
                            on.orf = TRUE,
                            ignore.illegal = TRUE
                            )

                  demo_DNA.seq <- seqinr::read.fasta(fastapath)

                  find_ORF <- getFromNamespace("find_ORF", "LncFinder")
                  Internal.EucDistance <- getFromNamespace("Internal.EucDistance", "LncFinder")
                  Internal.LogDistance <- getFromNamespace("Internal.LogDistance", "LncFinder")
                  Internal.hexamerScore <- getFromNamespace("Internal.hexamerScore", "LncFinder")

                  seq_orf <- lapply(demo_DNA.seq, function(x) {
                    orf.info <- find_ORF(x, max.only = TRUE)
                    if (orf.info[[1]] >= 12) seqinr::s2c(orf.info[[3]]) else NA
                  })
                  seq_orf <- seq_orf[!is.na(seq_orf)]

                  EucDis <- data.frame(t(sapply(
                    seq_orf,
                    Internal.EucDistance,
                    k = 6,
                    step = 1,
                    alphabet = c("a", "c", "g", "t"),
                    referFreq = referFreq
                  )))
                  LogDistance <- data.frame(t(sapply(
                    seq_orf,
                    Internal.LogDistance,
                    k = 6,
                    step = 1,
                    alphabet = c("a", "c", "g", "t"),
                    referFreq = referFreq
                  )))
                  hexamerScore <- data.frame(Hexamer.Score = sapply(
                    seq_orf,
                    Internal.hexamerScore,
                    k = 6,
                    step = 1,
                    alphabet = c("a", "c", "g", "t"),
                    referFreq = referFreq
                  ))

                hdata2<-cbind(EucDis,LogDistance)
                result <- cbind(hdata2,hexamerScore)
                return(result)

                }'''
              )

    sstruc = robjects.r['makeORFEucDist'](fastapath, abs_cds_path, abs_lncrna_path)
    return sstruc
