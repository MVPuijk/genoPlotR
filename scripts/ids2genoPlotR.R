#!/usr/bin/env Rscript

library(genoPlotR)
library(ape)
library(parallel)
cat("\n")


usage <- '

  USAGE:

      Example 1 (basic):
              Rscript ids2genoplotR.R gbk_files=*gbk ids=ids_file.tsv

      Example 2 (basic + comparisons + tree):
              Rscript ids2genoplotR.R gbk_files=*gbk comp_files=*comp comp_format=tab ids=ids_file.tsv tree=tree.nwk

      Example 3 (basic + comparisons + tree + annotations + longer (20 kb) neighbourhood + print positions to file + large PDF output + output file name):
              Rscript ids2genoplotR.R gbk_files=*gbk comp_files=*comp comp_format=tab ids=ids_file.tsv tree=tree.nwk annots=yes print_positions=yes xlim_out=positions.xlim out=genoplotr.pdf length=20000 heigth=20 width=30


  Arguments:

      - ids (Mandatory): tab-delimited file with the IDs of the genes to be printed, and the following structure: \'ID Genome_name color plot\'
         - ID: value of the genbank qualifyers \'synonym\' OR \'proteinid\'
         - Genome_name: genome name where to look for the ID (affects search of gbk_files and comp_files)
         - Color: anything identifiable as a color by R (e.g. "red", "blue4", "salmon")
         - Plot: optional field. If the term "plot" is present, ids2genoplotR will plot the genome neighbourhood of the given ID. If absent, the script will simply color the gene if it happens to be present in the neighbourhood delimited by a different ID

      - gbk_files (Mandatory): genbank files. Multi-contig genomes need to be concatenated in a contiguous genbank file, and contigs need to be encoded as features with the key \'contig\'


      - length: neighbourhood size to plot, unless a contig boundary found earlier (bp) (default: 10000).

      - comp_files: comparison files for the genomes used. The file names should contain the names of the genomes joined by an underscore, e.g. \'Genome1_Genome2.suffix\'

      - comp_format: choose between \'blast\' or \'tab\'
                - \'blast\': blast comparison file (blastN or tblastX)
                - \'tab\': a custom tab-delimited file with the following structure: \'start1 end1 start2 end2 col gene1 gene2\'

      - out: name of output pdf file (default: \'out.pdf\')

      - tree: newick file including as leaves the genome names used above

      - height: pdf height, in inches (default: 7)

      - width: pdf width, in inches (default: 20)

      - annotations/annots: if \'yes\' or \'1\', gene anontation (\'gene name\' and \'product\') will be plotted

      - print_positions: if \'yes\' or \'1\', neighbourhood limit positions will be printed to screen

      - xlim_out: name of the xlim output file, if "print_positions" is used. If not specified, xlim will be printed to a file name as the output file, with suffix ".xlim".  

      - xlim: input file with predefined positions to print, as outputted by "print_positions". One line per genome, in plotting order, with positions separated by tabs, e.g.: \'1 10000 155000 154000\' 

      - offset: input file with predefined offset values to plot

      - threads/nt: number of threads for parallel processing (default: 1)


'

#####################################################################
#####################################################################
#                                                                   #
# Print usage if no arguments are given. Parse arguments otherwise  #
#                                                                   #
#####################################################################
#####################################################################

# Functions for printing to screen 
col_start <- "\033[0;0;32m"
col_end <- "\033[0m"
catgreen <- function(text) {
	cat(paste0("\033[0;32m", text, col_end, "\n"))
}
catyellow <- function(text) {
	cat(paste0("\033[0;33m", text, col_end, "\n"))
}

# Read arguments

# Test arguments
#a <- "gbk_files=genomes_genoplotr/*concat comp_files=blast/*comp.2 comp_format=tab ids=ids_T3SS.col.noGa.test3 out=z_test.pdf tree=../pvc_selected_lineages.fasta.treefile.ed.ordered.renamed print_positions=yes height=20"
#args <- unlist(strsplit(a," "))

args <- commandArgs(TRUE)
if (length(args) == 0) {
	options(warning.length = 8170)
	stop(usage)
}

# Check if any provided arguments are unreadable
do_not_match <- grep("^gbk_files=|^ids=|^length=|^comp_files=|^comparison_files=|^comp_format=|^comparison_format=|^out=|^annots=|^annotations=|^tree=|^height=|^width=|^print_positions=|^print_pos=|^xlim=|^offset=", args, invert=TRUE)
if (length(do_not_match) > 0) {
	catyellow("#####################################################\n")
	catyellow(    paste("\nCould not understand the following arguments:\n", paste(args[do_not_match],collapse="\n"), "\n\nPlease check your command.\nYou can run this script without arguments to read the usage information:\n Rscript /local/one/dtamarit/scripts/ids2genoPlotR.R\n\n", sep=""))
	catyellow("#####################################################\n")
	stop()
}

# Read gbk files
gbk_files_arg <- grep("^gbk_files=", args)
gbk_files <- unlist(strsplit( args[ gbk_files_arg ], "="))[2]
filenames = Sys.glob(gbk_files)
if (! length(filenames)) {
	stop("Did not find gbk files")
}

# Read ID file, and prepare the list of genomes to print
ids_arg <- grep("^ids=", args)
ids_file <- unlist(strsplit( args[ ids_arg ], "="))[2]
ids <- read.csv( ids_file, header=F, col.names=c("id","genome","col","type"), stringsAsFactors=F, sep="\t" )
ids$genome = as.factor(ids$genome)
genomes_print = sort(unique(ids[ids$type=="plot","genome"]))

filenames = sort(filenames[grep(paste(genomes_print,collapse="|"),filenames)])
ids = ids[grep(paste(genomes_print,collapse="|"),ids[,"genome"]),]

# Fail-safe: If a genome is provided in the ID file but its genome was not found, exit with an error
for (i in 1:length(genomes_print)) {
	if (! length(grep(genomes_print[i],filenames))) {
		stop(paste("\n",genomes_print[i], " (found in IDs file) was not found in filenames (gbk_files):\n", paste(filenames, collapse="\n"), "\n\n", sep=""))
	}
}

# Read comparison files argument
use_comparisons = 0
comparison_arg <- grep("^comparison_files=|^comp_files=", args)
comparison_files <- unlist(strsplit( args[ comparison_arg ], "="))[2]
if (length(comparison_files)) {
	use_comparisons = 1
	comparison_filenames=Sys.glob(comparison_files)
	comp_format_arg <- grep("^comparison_format|^comp_format=", args)
	if (length(comp_format_arg)) {
    	comp_format <- unlist(strsplit( args[ comp_format_arg ], "="))[2]
		if ( (comp_format != "blast") & (comp_format != "tab") ) {
			stop(paste0("Could not understand comp_format argument as <",comp_format,">.\n\nPlease check your command.\nYou can run this script without arguments to read the usage information:\n Rscript /local/one/dtamarit/scripts/ids2genoPlotR.R\n\n"))
		}
	} else {
		comp_format = "blast"
	}
	catgreen(paste0("Will read comparison files as ", comp_format, "\n"))
}

# Read output PDF argument
out_arg <- grep("^out=", args)
outfile <- unlist(strsplit( args[ out_arg ], "="))[2]
if (length(outfile)) {
	if (file.exists(outfile)) {
		catgreen(paste0("Output file ",outfile," exists. Overwriting...\n") )
	}
}
if (length(outfile) == 0) {
	outfile="out.pdf"
}
catgreen(paste("Plotting to ",outfile,"\n",sep=""))


# Read whether to use annotations
annots_arg <- grep("^annots=|^annotations=", args)
annots_param <- unlist(strsplit( args[ annots_arg ], "="))[2]
if(length(annots_param)) {
	if ( annots_param=="yes" || annots_param==1 ) {
		use_annots=1
	} else {
		use_annots=0
	}
} else {
	use_annots=0
}

# Read tree and trim unused leaves
tree_arg <- grep("^tree=", args)
treefile <- unlist(strsplit( args[ tree_arg ], "="))[2]
if (length(treefile)) {
	if ( file.exists(treefile) ) {
		catgreen(paste0("Reading tree from ",treefile, "\n"))
		tree_read=read.tree( treefile )
		if ( length(tree_read$tip.label) > length(genomes_print) ) {
			remove = tree_read$tip.label[grep(paste(as.character(genomes_print),collapse="|"),tree_read$tip.label,invert=TRUE)]
			tree_read =  drop.tip(tree_read,remove)
		}
		write.tree(tree_read,"___outtree.tmp")
		newick <- readLines("___outtree.tmp")
		file.remove("___outtree.tmp")
		tree = newick2phylog(newick)
		dna_seg_order = tree_read$tip.label
		ids$genome = factor(ids$genome, levels=dna_seg_order)
		ids = ids[ order(ids$genome),]
		genomes_print = sort(factor(genomes_print, levels=dna_seg_order))
	} else {
		catgreen(paste0("Could not find treefile at ", treefile, ". Skipping tree...\n"))
		dna_seg_order = sort(genomes_print)
		ids = ids[ order(ids$genome),]
	}
	if (length(tree_read$tip.label) < length(genomes_print)) {
		stop("The given tree contains fewer taxa than required by the ID file. Check your data\n")
	}
} else {
	catgreen("No tree given; plotting without a tree...\n")
	dna_seg_order = sort(genomes_print)
	ids = ids[ order(ids$genome),]
}
ids$genome = as.character(ids$genome)
genomes_print = as.character(genomes_print)
dna_seg_order = as.character(dna_seg_order)

# Read length of printed genome regions. 10 kb by default 
length_arg <- grep("^length=", args)
if ( length(length_arg) ) {
	if ( grep("^[0-9]+$", length_arg, perl=TRUE) ) {
		length <- as.numeric(unlist(strsplit( args[ length_arg ], "="))[2])
	} else {
		stop(paste0("Argument 'length' (",length,") was not recognized as numeric for segment length"))
	}
} else {
	length=10000
}
catgreen(paste("Plotting neighbourhoods up to ",length," bp from selected loci.\n",sep=""))

# Read plotting dimension parameters
height_arg <- grep("^height=", args)
height <- as.numeric(unlist(strsplit( args[ height_arg ], "="))[2])
if (length(height) == 0) {
#  height = 7
   height = "calculate"
   catgreen("No plot height argument was provided. It will be calculated based on the plotted data\n")
}

width_arg <- grep("^width=", args)
width <- as.numeric(unlist(strsplit( args[ width_arg ], "="))[2])
if (length(width) == 0) {
#  width = 20
   width = "calculate"
   catgreen("No plot width argument was provided. It will be calculated based on the plotted data\n")
}

# Read xlim outfile name 
xlim_out_arg <- grep("^xlim_out=", args)
xlim_out_param <- unlist(strsplit( args[ xlim_out_arg ], "="))[2]
xlimfile=sub(".pdf", ".xlim",outfile)
if (length(xlim_out_param)) {
	xlimfile=xlim_out_param
} 

# Read whether to print positions
print_pos_arg <- grep("^print_positions=|^print_pos=", args)
print_pos_param <- unlist(strsplit( args[ print_pos_arg ], "="))[2]
print_positions = 0
if( length(print_pos_param) ) {
	if ( print_pos_param=="yes" || print_pos_param==1 ) {
		print_positions = 1
		catgreen(paste0("Printed positions will be printed to screen and to file ",xlimfile, "\n"))
	}
}

# Read xlim input file
xlim_arg <- grep("^xlim=",args)
xlim_param <- unlist(strsplit( args[ xlim_arg ], "="))[2]
if (length(xlim_param)) {
	if (file.exists(xlim_param)) {
		catgreen(paste0("xlim file detected. Will plot positions as defined by ",xlim_param,"\n"))
	} else {
		stop(paste0("Could not find xlim file ",xlim_param, "\n"))
	}
}

# Read offset input file
offset_arg <- grep("^offset=",args)
offset_param <- unlist(strsplit( args[ offset_arg ], "="))[2]
if (length(offset_param)) {
	if (file.exists(offset_param)) {
		catgreen(paste0("offset file detected. Will plot positions as defined by ",offset_param,"\n"))
	} else {
		stop(paste0("Could not find offset file ",xlim_param, "\n"))
	}
}

# Read number of threads -- CURRENTLY NOT IMPLEMENTED
threads_arg <- grep("^threads=|^nt=", args)
threads <- as.numeric(unlist(strsplit( args[ threads_arg ], "="))[2])
if (length(threads) == 0) {
	threads = 1
}

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################


#####################################################################
#####################################################################
#                                                                   #
#        Read genomes and edit basic plotting parameters            #
#                                                                   #
#####################################################################
#####################################################################



catyellow(paste("\nRunning genoPlotR with ", length(genomes_print), " genomes:\n", sep=""))
cat(paste(genomes_print,sep=","))
cat("\n\n")

catyellow("Reading genomes...\n")
for (i in 1:length(genomes_print)) {
	file = filenames[grep(genomes_print[i], filenames)]
	cat(paste("   Reading ", file, "\n", sep=""))

	eval(bquote(
		.(as.name(genomes_print[i])) <- read_dna_seg_from_genbank(file, tagsToParse = c("CDS", "tRNA", "rRNA","repeat_region","contig"))
	))
}

catyellow("\nProcessing genomes...\n")
for (i in 1:length(genomes_print)) {
	eval(bquote(  .(as.name(genomes_print[i]))[["col"]] <- "black"  ))
	eval(bquote(  .(as.name(genomes_print[i]))[["fill"]] <- "gray95"  ))
	eval(bquote(  .(as.name(genomes_print[i]))[["gene_type"]] <- "lines"  ))
	eval(bquote(  .(as.name(genomes_print[i]))[ .(as.name(genomes_print[i]))$feature == "CDS", "gene_type" ] <- "arrows"  ))
	eval(bquote(  .(as.name(genomes_print[i]))[ .(as.name(genomes_print[i]))$feature == "tRNA", "gene_type" ] <- "blocks"  ))
	eval(bquote(  .(as.name(genomes_print[i]))[ .(as.name(genomes_print[i]))$feature == "rRNA", "gene_type" ] <- "headless_arrows"  ))
	eval(bquote(  .(as.name(genomes_print[i]))[ .(as.name(genomes_print[i]))$feature == "repeat_region", "gene_type" ] <- "blocks"  ))
	eval(bquote(  .(as.name(genomes_print[i]))[ .(as.name(genomes_print[i]))$feature == "tRNA", "fill" ] <- "darkseagreen1"  ))
	eval(bquote(  .(as.name(genomes_print[i]))[ .(as.name(genomes_print[i]))$feature == "rRNA", "fill" ] <- "lightblue1"  ))
	eval(bquote(  .(as.name(genomes_print[i]))[ .(as.name(genomes_print[i]))$feature == "repeat_region", "fill" ] <- "lemonchiffon"  ))
}





#####################################################################
#####################################################################
#                                                                   #
#    Loop over IDs,        	                                        #
#    modify their colors,         	                                #
#    obtain positions of their neighbourhoods,        	            #
#    and create objects to denote contig boundaries                 #
#                                                                   #
#####################################################################
#####################################################################


catyellow("\nObtaining neighbourhoods from gene IDs...\n")
limits = data.frame("genome" = 1, "start" = 1, "end" = 1)
for (i in 1:length(ids[,"id"])) {
	# Change colors for each provided ID
	eval(bquote(
		.(as.name(  ids[i,"genome"]  ))[
			grepl( paste0(  ids[i,"id"]  ,"$"), .( as.name(  ids[i, "genome"]  ))[['synonym']]) |
			grepl( paste0(  ids[i,"id"]  ,"$"), .( as.name(  ids[i, "genome"]  ))[['proteinid']])
			,
			 "fill" ] <- ids[i,"col"]
	))
	if ( ids[i, "type"] != "plot") next

	# Get positions for each ID whose neighbourhoods are to be plotted  
	line = eval(bquote( .(as.name(  ids[i,"genome"]  ))[
			grepl( paste(  ids[i,"id"]  ,"$",sep=""), .( as.name(  ids[i, "genome"]  ))[['synonym']]) |
			grepl( paste(  ids[i,"id"]  ,"$",sep=""), .( as.name(  ids[i, "genome"]  ))[['proteinid']])
			,
			]
		))

	if (length(line$name) == 0 ) {
		cat(paste("Could not find ID ", ids[i,"id"], " in genome ", ids[i, "genome"], "\n", sep = ""))
		next
	}
	limits[i, "genome"] = ids[i, "genome"]
	limits[i, "start"] = line$start
	limits[i, "end"] = line$end
}
# Clean up and order appropriately
limits = limits[ complete.cases(limits),]
limits = limits[ limits$genome != "1", ]
limits2 = limits[0,]
for (g in dna_seg_order) {
    grepped = limits[grep(g, limits$genome),]
    grepped = grepped[ with (grepped, order(start)), ]
    limits2 = rbind(limits2, grepped)
}
limits = limits2
catyellow("\nNeighbourhood limits calculated;")


# Generate xlim object and generate features representing contig boundaries 
xlim = list()
for (i in 1:length(limits[,1])) {

	# Obtain contig start and end positions
	
	start <- limits[i,"start"]
	end <- limits[i,"end"]

	contig_start = eval(bquote(
		.(as.name( limits[i,"genome"] ))[ which(
			.(as.name( limits[i,"genome"] ))$feature == "contig" &
			.(as.name( limits[i,"genome"] ))$start <= start &
			.(as.name( limits[i,"genome"] ))$end >= end
			)
			, "start" ]
	))

	contig_end = eval(bquote(
		.(as.name( limits[i,"genome"] ))[ which(
			.(as.name( limits[i,"genome"] ))$feature == "contig" &
			.(as.name( limits[i,"genome"] ))$start <= start &
			.(as.name( limits[i,"genome"] ))$end >= end
			)
			, "end" ]
	))

	# Determine gene neighbourhood limits per gene of interest: either a fixed distance from gene start, or contig boundary
	# Add feature lines to genome objects representing contig boundaries
	
	if (length(contig_start) == 1) {
		xlim_start = max( (start - length),   contig_start )
		if (xlim_start == contig_start) {
			eval(bquote(
				.(as.name( limits[i,"genome"] ))  <-
					rbind(
						.(as.name( limits[i,"genome"] ))
						,
						c("NA", contig_start+1, contig_start+9, 1, 3, "NA", "-", "NA", "NA", "NA", "CDS", "bars", "black", "black", 1, 3, 8, 1)
            	     )
			))
		}
	} else {
		xlim_start = (start - length)
	}
	
	if (length(contig_end) == 1) {
		xlim_end =   min( ( end  + length),   contig_end )
		if (xlim_end == contig_end) {
			eval(bquote(
			.(as.name( limits[i,"genome"] ))  <-
				rbind(
					.(as.name( limits[i,"genome"] ))
					,
					c("NA", contig_end-9, contig_end-1, 1, 3, "NA", "-", "NA", "NA", "NA", "CDS", "bars", "black", "black", 1, 3, 8, 1)
				)
			))
		}
	} else {
		xlim_end =  ( end  + length)
	}

	# Format
	eval(bquote(
		.(as.name( limits[i,"genome"] ))$start <- as.numeric( .(as.name( limits[i,"genome"] ))$start  )
	))
	eval(bquote(
		.(as.name( limits[i,"genome"] ))$end <- as.numeric( .(as.name( limits[i,"genome"] ))$end  )
	))
	eval(bquote(
		.(as.name( limits[i,"genome"] ))$lwd <- as.numeric( .(as.name( limits[i,"genome"] ))$lwd  )
	))
	eval(bquote(
		.(as.name( limits[i,"genome"] ))$lty <- as.numeric( .(as.name( limits[i,"genome"] ))$lty  )
	))
	eval(bquote(
		.(as.name( limits[i,"genome"] ))$pch <- as.numeric( .(as.name( limits[i,"genome"] ))$pch  )
	))
	eval(bquote(
		.(as.name( limits[i,"genome"] ))$cex <- as.numeric( .(as.name( limits[i,"genome"] ))$cex  )
	))


	# Extend a previous gene neighbourhood if overlapping
	genome_number = grep( limits[i,"genome"],genomes_print )
	if (i > 1) {
		if (  genomes_print[genome_number]  !=  limits[i-1,"genome"] ) {

			xlim[[ genome_number  ]] <- c( xlim_start, xlim_end )

		} else {
			if (  xlim[[ genome_number ]][length(xlim[[ genome_number ]])] > xlim_start ) {
				xlim[[ genome_number ]][length(xlim[[ genome_number ]])] <- xlim_end
			} else {
				xlim[[ genome_number ]] <- c(xlim[[ genome_number ]], xlim_start, xlim_end )
			}
		}
	} else {
		xlim[[ genome_number  ]] <- c( xlim_start, xlim_end )
	}

}
catyellow(" contig boundaries formatted for plotting.\n")





#####################################################################
#####################################################################
#                                                                   #
#           Detect if positions to print were provided              #
#           and overwrite xlim parameter                            #
#                                                                   #
#####################################################################
#####################################################################

if (length(xlim_param)) {

	xlim_input <- list()
	xlim_read <- read.delim(xlim_param,sep=" ",header=F)

	old_i = 0
	for(i in 1:nrow(xlim_read)) {
		for(j in 1:ncol(xlim_read)) {
			if (! is.na(as.numeric(xlim_read[i,j]))) {
				if(i != old_i) {
					xlim_input[[ i ]] <- xlim_read[i,j]
				} else {
					xlim_input[[ i ]] <- c(xlim_input[[ i ]], xlim_read[i,j])
				}
			}
			old_i <- i
		}
   }
   xlim <- xlim_input
   catyellow(paste0("xlim parameter obtained from input file ",xlim_param,".\n"))
}

catyellow("xlim parameter ready\n\n")



#####################################################################
#####################################################################
#                                                                   #
#         Print obtained neighbourhood positions if asked           #
#                                                                   #
#####################################################################
#####################################################################


if (print_positions) {
	catyellow("Printing neighbourhood limits:\n\n")
	cat("--------------------------------\n")
	for (i in 1:length(genomes_print)) {
		cat(genomes_print[i])
		cat("\n")
		for ( j in 1:( length(xlim[[i]])/2 ) ) {
			cat( xlim[[i]][j*2-1] )
			cat(" ")
			cat( xlim[[i]][j*2] )
			cat("\n")
		}
		cat("\n")
	}
	cat("--------------------------------\n")
	cat("\n")

#  xlimfile <- sub(".pdf",".xlim",outfile)
	if (file.exists(xlimfile)) {
		file.remove(xlimfile)
	}
	invisible(lapply(xlim,write,xlimfile,append=TRUE,ncolumns=1000))
  
}


#####################################################################
#####################################################################
#                                                                   #
#              Generate dna_segs and prepare annotations            #
#                                                                   #
#####################################################################
#####################################################################

dna_segs = list()
annot = list()
for (i in 1:length(genomes_print)) {

	eval(bquote(  dna_segs[[i]] <- as.dna_seg( .(as.name(genomes_print[i])) ) ))

	annot[[i]] <- annotation(
                  x1=( dna_segs[[i]][dna_segs[[i]]$feature == "CDS", "start"] + dna_segs[[i]][dna_segs[[i]]$feature == "CDS", "end"] ) / 2,
                  rot=ifelse(i==1,-30,30),
                  text=
                       paste(
                          dna_segs[[i]][dna_segs[[i]]$feature == "CDS", "synonym"],
#                          dna_segs[[i]][dna_segs[[i]]$feature == "CDS", "proteinid"],
                          sub("hypothetical protein", "", dna_segs[[i]][dna_segs[[i]]$feature == "CDS", "product"])
                        , sep="  ")
            )
}

#####################################################################
#####################################################################
#                                                                   #
#                Read comparison files if asked                     #
#                                                                   #
#####################################################################
#####################################################################

if (use_comparisons) {
	catyellow("Reading comparison files...\n")

	comparisons = list()
#  if (length(treefile)) {
#    genome_print = sort(tree_read$tip.label)
#  }
	for (i in 1:(length(genomes_print)-1)) {
		comp_file = comparison_filenames[grep(paste(genomes_print[i],"_",genomes_print[i+1],sep=""), comparison_filenames)]
		cat(paste("   Reading ",comp_file, "\n", sep=""))
		if (comp_format == "blast") {
			comparisons[[i]] = read_comparison_from_blast( comp_file )
		}
		if (comp_format == "tab") {
			comparisons[[i]] = read_comparison_from_tab( comp_file ,header=TRUE)
		}
	}
	cat("\n")
} else  {
	catyellow("No comparison files specified. Plotting without comparisons...\n\n")
}

#####################################################################
#####################################################################
#                                                                   #
#           Detect if positions for offsets parameter               #
#           were provided and parse them                            #
#                                                                   #
#####################################################################
#####################################################################

if (length(offset_param)) {

	offset_input <- list()
	offset <- as.numeric(readLines(offset_param))

   catyellow(paste("offsets parameter obtained from input file ",offset_param,".\n\n\n",sep=""))
} else {
	offset <- ""
}






##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#####################################################################
#####################################################################
#                                                                   #
#           Calculate plotting parameters, if undefined             #
#                                                                   #
#                                                                   #
#####################################################################
#####################################################################

# Calculate plot height as a function of the number of genomes and whether comparisons or annotations are plotted
if (height == "calculate") {
	height = length(dna_segs)/2
	if (use_annots) {
		height = height*4
	} else {
		if (use_comparisons) {
			height = height*2
		}
	}
	catyellow(paste0("Calculated plot height: ", height, "\n"))
}

# Calculate plot width as a function of the maximum number of bases to be printed for a genome
if (width == "calculate") {
	len = numeric(length(xlim))
	for (i in 1:length(xlim)) {       
		for (j in 1:(length(xlim[[i]])/2)) {
			start <- xlim[[i]][j*2-1]
			end <- xlim[[i]][j*2]
			len[i] <- len[i] + max(start,end)-min(start,end)
		}
	}
	width = max(len)
	width = width / 4000

	if (length(treefile)) {
		width = width + 5
	}

	if (width < 5) {
		width = 5
	}

	catyellow(paste0("Calculated plot width: ", width, "\n"))
}





#####################################################################
#####################################################################
#                                                                   #
#                            Plot                                   #
#                                                                   #
#####################################################################
#####################################################################


catyellow(paste("Plotting (h=",height,", w=",width,") at ",outfile, "\n",sep=""))

if ((! use_comparisons) & (! length(treefile)) & (! use_annots) ) {
	catyellow("Only plotting dna_segs\n\n")
	pdf(outfile, height=max(height,length(genomes_print)/height), width=width)
	plot_gene_map(
		dna_segs = dna_segs,
		dna_seg_labels = dna_seg_order,
		xlims = xlim
	)
	dev.off()
}

if (use_comparisons & (! length(treefile)) & (! use_annots)) {
	catyellow("Including:\n-comparisons\n\n")
	pdf(outfile, height=max(height,length(genomes_print)/height), width=width)
	plot_gene_map(
		dna_segs = dna_segs,
		dna_seg_labels = dna_seg_order,
		comparisons = comparisons,
		xlims = xlim
	)
	dev.off()
}

if ((! use_comparisons) & length(treefile) & (! use_annots)) {
	catyellow("Including:\n-tree\n\n")
	pdf(outfile, height=max(height,length(genomes_print)/height), width=width)
	plot_gene_map(
		dna_segs = dna_segs,
		dna_seg_labels = dna_seg_order,
		tree=tree,
		xlims = xlim
	)
	dev.off()
}

if ((! use_comparisons) & (! length(treefile)) & use_annots) {
	catyellow("Including:\n-annotations\n\n")
	pdf(outfile, height=max(height,length(genomes_print)/height), width=width)
	plot_gene_map(
		dna_segs = dna_segs,
		dna_seg_labels = dna_seg_order,
		annotations = annot,
		xlims = xlim
	)
	dev.off()
}


if (use_comparisons & length(treefile) & (! use_annots)) {
	catyellow("Including:\n-comparisons\n-tree\n\n")
	pdf(outfile, height=max(height,length(genomes_print)/height), width=width)
	plot_gene_map(
		dna_segs = dna_segs,
		dna_seg_labels = dna_seg_order,
		tree=tree,
		comparisons = comparisons,
		xlims = xlim
	)
	dev.off()
}

if ((! use_comparisons) & length(treefile) &  use_annots) {
	catyellow("Including:\n-tree\n-annotations\n\n")
	pdf(outfile, height=max(height,length(genomes_print)/height), width=width)
	plot_gene_map(
		dna_segs = dna_segs,
		dna_seg_labels = dna_seg_order,
		tree=tree,
		annotations = annot,
		xlims = xlim
	)
	dev.off()
}
if (use_comparisons & (! length(treefile)) & use_annots) {
	catyellow("Including:\n-comparisons\n-annotations\n\n")
	pdf(outfile, height=max(height,length(genomes_print)/height), width=width)
	plot_gene_map(
		dna_segs = dna_segs,
		dna_seg_labels = dna_seg_order,
		comparisons = comparisons,
		annotations = annot,
		xlims = xlim
	)
	dev.off()
}

if (use_comparisons & length(treefile) & use_annots) {
	catyellow("Including:\n-comparisons\n-tree\n-annotations\n\n")
	pdf(outfile, height=max(height,length(genomes_print)/height), width=width)
	plot_gene_map(
		dna_segs = dna_segs,
		dna_seg_labels = dna_seg_order,
		tree=tree,
		comparisons = comparisons,
		annotations = annot,
		xlims = xlim
	)
	dev.off()
}

###### TO BE FIXED: ONLY CURRENT IMPLEMENTATION OF OFFSETS PARAMETER

if (use_comparisons & (! length(treefile)) & (! use_annots) & ( length(offset_param)) ) {
	cat("Including:\n-comparisons\n-offsets\n\n")
	pdf(outfile, height=max(height,length(genomes_print)/height), width=width)
	plot_gene_map(
		dna_segs = dna_segs,
		dna_seg_labels = dna_seg_order,
		comparisons = comparisons,
		offsets=offset,
		xlims = xlim
	)
	dev.off()
}