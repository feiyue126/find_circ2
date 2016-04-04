#!/usr/bin/env Rscript

#######################
# general information #
#######################

# file:         generate_synthetic_reads.r
# created:      2016-03-21
# last update:  2016-04-04
# author:       Marcel Schilling <marcel.schilling@mdc-berlin.de>
# purpose:      generate synthetic reads based on known circRNAs (in BED format)


######################################
# change log (reverse chronological) #
######################################

# 2016-04-04: fixed read naming
#             fixed collapsing of all fragments from the same circRNA to the same start position /
#             removed duplicates
# 2016-03-29: switched to relative coordinates for splice junction annotation as requested by Marvin
#             added merging of mate information as requested by Marvin
# 2016-03-21: initial version (uncommented messy script with fixed parameters for EA_cel10T sample)


#############
# libraries #
#############

#require(doMC) # parallelize plyr calls
require(AnnotationHub)
require(magrittr)
require(GenomicRanges)
require(plyr)
require(GenomicFeatures)


##############
# parameters #
##############

#threads<-detectCores()
bsgenome<-"BSgenome.Celegans.Ensembl.WBcel235"
annotation<-"Caenorhabditis_elegans.WBcel235.83.gtf"
circs.bed<-"wbcel235.circs.bed.gz"
contigs.to.keep<-c("I","II","III","IV","V","X")
colnames.circs.bed<-c("chr"
                     ,"start"
                     ,"end"
                     ,"name"
                     ,"n_reads"
                     ,"strand"
                     ,"n_uniq"
                     ,"uniq_bridges"
                     ,"best_qual_left"
                     ,"best_qual_right"
                     ,"tissues"
                     ,"tiss_counts"
                     ,"edits"
                     ,"anchor_overlap"
                     ,"breakpoints"
                     ,"signal"
                     ,"strandmatch"
                     ,"category"
                     )
nreads_min<-10
sep.coord<-","
sep.fasta<-"|"
sep.isform<-":"
reads.circ<-1000
fragment.length<-300
sep.annotation<-"___"
read.length<-100
tag.origin<-"O"
tag.match<-"M"
tag.linear_splice<-"LS"
tag.circular_splice<-"CS"
sep.tag<-":"
sep.tags<-";"
sep.mates<-"|"
read1.fa<-"synthetic_reads.R1.fa.gz"
read2.fa<-"synthetic_reads.R2.fa.gz"


##################
# initialization #
##################

#registerDoMC(threads) # initialize parallelization engine for plyr


###############
# load genome #
###############

require(bsgenome,character.only=T)

# alias genome object
eval(parse(text=paste("bsgenome <-",bsgenome)))


###################
# load annotation #
###################

# load annotation hub GRanges
annotations<-AnnotationHub() %>%
             query("GRanges")

# try querying annotation to use by name
if(annotation %in% names(annotations)){
  annotation<-annotations[[annotation]]

# in case of failure, try querying by title
} else {
  annotations<-annotations %>%
               subset(title==annotation)
      
  # if successful, get annotation GRanges object
  if(length(annotations)){
    annotation<-annotations[[1]]

  # if unsuccessful, exit with error
  } else {
    annotation %>%
    paste0("error: Could not find GRanges object for annotation '"
          ,.
          ,"'."
          ) %>%
    stop
    
  }
}

# clean up annotation hub
rm(annotations)

seqlevels(annotation,force=T)<-contigs.to.keep


#################
# load circRNAs #
#################

# read find_circ output & filter out circRNAs
circs<-circs.bed %>%
       read.table %>%
       setNames(colnames.circs.bed) %>%
       {.[grep("CIRCULAR",.$category),]} %>%
       {.[.$chr %in% contigs.to.keep,]} %$%
       {.[.$n_reads>=nreads_min,]} %$%
       {
         GRanges(chr
                ,IRanges(start+1,end)
                ,strand) %>%
         setNames(name)
       }


#########################
# get circRNA sequences #
#########################

circs %<>%
llply(function(circ)
        circ %>%
        subsetByOverlaps(annotation,.) %>%
        makeTxDbFromGRanges %>%
        (
          function(txdb){
            txdb %>%
            transcripts %>%
            {subset(.,start(.)<=start(circ) & end(.)>=end(circ))} %$%
            {
              txdb %>%
              exonsBy('tx',use.names=T) %>%
              `[`(tx_name)
            }
          }
        ) %>%
        llply(. %>%
              {
                start(.)[which.min(start(.))]<-start(circ)
                end(.)[which.max(end(.))]<-end(circ)
                .
              }
#             ,.parallel=T
             ) %>%
        GRangesList %>%
        {
          if(length(.)){
            .
          } else {
            circ %>% 
            GRangesList
          }
        }
#     ,.parallel=T
     ) %>%
llply(function(isoforms)
        isoforms %>%
        extractTranscriptSeqs(bsgenome,.) %>%
        setNames(isoforms %>%
                 sort %>%
                 llply(.%$%
                       {
                         seqnames %>%
                         paste(start %>%
                               `-`(1) %>% 
                               paste0(collapse=sep.coord)
                              ,end %>%
                               paste0(collapse=sep.coord)
                              ,strand
                              ,sep=sep.fasta
                              ) %>%
                         unique
                       }
                      ) %>%
                 unlist
                )
     ) %>%
llply(unique) %>%
(
  function(circs)
    circs %>%
    names %>%
    llply(function(circ){
            circ %>%
            circs[[.]] %>%
            {
              if(length(.)>1){
                1:length(.) %>%
                paste(circ
                     ,.
                     ,sep=sep.isform
                     )
              } else {
                circ
              }
            } %>%
            paste(circs[[circ]] %>%
                  names
                 ,sep=sep.fasta
                 ) %>%
            setNames(circs[[circ]],.)
          }
         )
) %>%
do.call(c,.)


######################
# generate fragments #
######################
fragments<-circs %>%
           sample(reads.circ
                 ,replace=T
                 ) %>%
           llply(function(circ)
                   circ %>%
                   length %>%
                   seq(from=1) %>%
                   sample(1) %>%
                   (
                     function(startpos){
                       while(length(circ)<(startpos+fragment.length))
                         circ %<>% append(circ)
                       circ %>%
                       subseq(startpos
                             ,width=fragment.length
                             ) %>%
                       list %>%
                       setNames(startpos)
                     }
                   )
#                ,.parallel=T
                ) %>%
           {
             setNames(.
                     ,llply(.,names) %>%
                      unlist %>%
                      paste(names(.)
                           ,.
                           ,sep=sep.annotation
                           )
                     )
           } %>%
           llply(unname) %>%
           unlist %>%
           llply(toString) %>% # This should not be necessary,
           unlist %>%          # but it is!
           DNAStringSet %>%
           unique




#############
# get reads #
#############

reads<-fragments %>%
       {
         name<-names(.) %>%
               strsplit(sep.annotation
                       ,fixed=T
                       )
         startpos<-name %>%
                   llply(`[`,2) %>%
                   unlist %>%
                   as.integer
         name %<>%
         llply(`[`,1) %>%
         unlist
         list(read1=subseq(.
                          ,1
                          ,read.length
                          )
             ,read2=subseq(.
                          ,fragment.length-read.length+1
                          ,fragment.length
                          ) %>%
                    setNames(name %>%
                             paste(startpos+fragment.length-read.length
                                  ,sep=sep.annotation
                                  )
                            )
             )
       }


##################
# annotate reads #
##################

reads %<>%
llply(function(reads){
        name<-reads %>%
              names %>%
              strsplit(sep.annotation
                      ,fixed=T
                      ) %>%
              do.call(rbind,.)
        startpos<-name %>%
                 `[`(,2) %>%
                  as.integer
        name %<>%
        `[`(,1) %>%
        strsplit(sep.fasta
                ,fixed=T
                ) %>%
        do.call(rbind,.)
        strand<-name[,5]
        ends<-name[,4] %>%
              strsplit(sep.coord
                      ,fixed=T
                      ) %>%
              llply(as.integer)
        starts<-name[,3] %>%
                strsplit(sep.coord
                        ,fixed=T
                        ) %>%
                llply(as.integer) %>%
                llply(`+`,1) # switch to 1-based closed intervals
        chr<-name[,2]
        name<-name[,1]
        starts %>%
        length %>%
        seq(1,.) %>%
        llply(function(read)
                list(chr=chr[read]
                    ,starts=starts[[read]]
                    ,ends=ends[[read]]
                    ,strand=strand[read]
                    ,name=name[read]
                    ,startpos=startpos[read]
                    ,read=reads[[read]] %>%
                          DNAStringSet
                    ,exons=starts[[read]] %>%
                           length
                    ,lengths=ends[[read]] %>%
                             `-`(starts[[read]]) %>%
                             `+`(1)
                    )
#             ,.parallel=T
             )
      }
     ) %>%
llply(llply
     ,. %$%
      {
        exon<-ifelse(strand=="+"
                    ,1
                    ,exons
                    )
        pregap<-0
        while((pregap+lengths[exon])<startpos){
          pregap %<>% `+`(lengths[exon])
          exon %<>% `+`(0-2*(strand=="-")) %>%
                    mod(exons) %>%
                    `+`(1)
        }
        readstarts<-ifelse(strand=="+"
                             ,starts[exon]+startpos-pregap-1
                             ,ends[exon]-startpos+pregap+1
                             )
        readends<-ifelse(strand=="+"
                        ,min(ends[exon],readstarts+read.length-1)
                        ,max(starts[exon],readstarts-read.length+1)
                        )
        readlen<-ifelse(strand=="+"
                       ,readends-readstarts+1
                       ,readstarts-readends+1
                       )
        readexon<-1
        while(readlen<read.length){
          exon %<>% `+`(0-2*(strand=="-")) %>%
                    mod(exons) %>%
                    `+`(1)
          readexon %<>% `+`(1)
          readstarts[readexon]<-ifelse(strand=="+"
                                      ,starts[exon]
                                      ,ends[exon]
                                      )
          readends[readexon]<-ifelse(strand=="+"
                                    ,min(ends[exon],readstarts[readexon]+read.length-readlen-1)
                                    ,max(starts[exon],readstarts[readexon]-read.length+readlen+1)
                                    )
          readlen %<>% `+`(ifelse(strand=="+"
                                 ,readends[readexon]-readstarts[readexon]+1
                                 ,readstarts[readexon]-readends[readexon]+1
                                 )
                          )
        }
        name %>%
        paste(tag.origin  %>%
              paste(chr
                   ,readstarts[1]-1 # switch back to 0-based coordinates
                   ,strand
                   ,sep=sep.tag
                   ) %>%
              paste(1:readexon %>%
                    llply(function(exon){
                            annot<-ifelse(strand=="+"
                                         ,readends[exon]-readstarts[exon]+1
                                         ,readstarts[exon]-readends[exon]+1
                                         ) %>%
                                   paste(tag.match
                                        ,.
                                        ,sep=sep.tag
                                        )
                            if(exon<readexon)
                              annot %<>% paste(ifelse(((strand=="+")
                                                      &(readends[exon]<readstarts[exon+1])
                                                      )
                                                     |((strand=="-")
                                                      &(readends[exon]>readstarts[exon+1])
                                                      )
                                                     ,tag.linear_splice
                                                     ,tag.circular_splice
                                                     ) %>%
                                               paste(readends[exon]-readstarts[1]     # switch to relative coordinates
                                                    ,readstarts[exon+1]-readstarts[1] # switch to relative coordinates
                                                    ,sep=sep.tag
                                                    )
                                              ,sep=sep.tags
                                              )
                            annot
                          }
                         ) %>%
                    unlist %>%
                    paste0(collapse=sep.tags)
                   ,sep=sep.tags
                   )
             ,sep=sep.annotation
             ) %>%
        setNames(read,.)
      }
#     ,.parallel=T
     ) %>%
llply(do.call
     ,what=c
#     ,.parallel=T
     ) %>%
(
  function(mates){
    names(mates$read1)<-mates %>%
                        llply(names) %$%
                        {

                          paste(read1
                               ,read2 %>%
                                strsplit(sep.annotation
                                        ,fixed=T
                                        ) %>%
                                `[[`(1) %>%
                                `[`(2)
                               ,sep=sep.mates
                               )
                        }
    names(mates$read2)<-names(mates$read1)
    mates
  }
)


###############
# write reads #
###############

writeXStringSet(reads$read1,read1.fa,compress=grepl("\\.gz$",read1.fa))
writeXStringSet(reads$read2,read2.fa,compress=grepl("\\.gz$",read1.fa))
