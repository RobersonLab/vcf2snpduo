#!/usr/bin/python

##########
#imports #
##########
import argparse
import sys
import time

########
# Fxns #
########
def log_msg( msg, printTime=True ):
	if printTime:
		sys.stderr.write( "%s\n" % ( time.strftime( '%H:%M %z on %b %d, %Y' ) ) )
	sys.stderr.write( "%s\n" % ( msg ) )
	sys.stderr.flush()
	
def error_msg( msg ):
	log_msg( "Error: %s\n" % ( msg.rstrip() ) )
	sys.exit( 1 )

def genotype_conversion( genotype ):
	if genotype == '0/0' or genotype == '0|0':
		return( "AA" )
	elif genotype == '0/1' or genotype == '0|1' or genotype == '1|0':
		return( "AB" )
	elif genotype == '1/1' or genotype == '1|1':
		return( "BB" )
	
	return( "NC" )

if __name__ == "__main__":
	####################
	# Version and name #
	####################
	SCRIPT_PATH = sys.argv[0]
	SCRIPT_NAME = SCRIPT_PATH.split( '/' )[-1].split( '\\' )[-1]
	VERSION = '1.0.2'
	
	#############
	# arg parse #
	#############
	parser = argparse.ArgumentParser()

	parser.add_argument( "input_vcf", help="VCF file to be processed" )
	parser.add_argument( "--minDepth", help="Minimum depth for individual samples", type=int, default=10, dest="minDepth" )
	parser.add_argument( "--minSampleCount", help="Minimum number of samples retained for a variant to print", type=int, default=1, dest="minCount" )
	parser.add_argument( "--filterDP", help="Default is filtering by individual Allele Depth. This option switches the filter to individual level DP. If neither are present, this tool isn't appropriate.", default=True, action='store_false', dest='filterByAD' )
	parser.add_argument( "--noDepthFilter", help="Default is allele depth (AD) filtering, with depth (DP) option. This option removes all depth filtering", default=False, action='store_true', dest='noFilter' )
	parser.add_argument( '--version', action='version', version="v%s" % ( VERSION ) )

	args = parser.parse_args()

	log_msg( "%s v%s\nMinimum Depth: %s\nMinimum Samples: %s\nDepth filter field: %s" % ( SCRIPT_NAME, VERSION, args.minDepth, args.minCount, "None" if args.noFilter else "AD" if args.filterByAD else "DP" ) )

	#####################
	# open file handles #
	#####################
	try:
		VCF = open( args.input_vcf , 'r' )
	except:
		error_msg( "Problem opening input [%s]" % ( args.input_vcf ) )

	colOrder = {}
	sampleNames = []
		
	# Prepare output header from sample names
	headerString = "Chromosome,Physical.Position,RSID,Ref,Var"
	for line in VCF:
		line = line.rstrip()
		
		if line[:6] == "#CHROM":
			vals = line.lstrip( '#' ).split( "\t" )
			
			formatSet = False
			
			for colIndex in range( len( vals ) ):
				colName = vals[colIndex]
				
				if colName in colOrder:
					error_msg( "Repeated header in VCF file" )
				else:
					colOrder[colName] = colIndex
					
					if colName == "FORMAT":
						formatSet = True
					elif formatSet == True:
						sampleNames.append( colName )
				
			if 'FORMAT' not in colOrder:
				error_msg( "Is this a sites file? Too few fields found." )
			else:
				print "%s,%s" % ( headerString, ','.join( sampleNames ) )
			
			break

	# process the actual data of the VCF
	for line in VCF:
		keepSampleCount = 0
		
		line = line.rstrip()
		lineVals = line.split( '\t' )
		#chrom, position, id, refBase, varBase, varScore, filter, detailInfo, sampleValueOrder, sampleData = line.split("\t", 9)
		
		# skip multi-allelic sites
		if len( lineVals[ colOrder['ALT'] ].split( "," ) ) > 1:
			log_msg( "ALT base [%s] has more than 1 value, skipping...\n" % ( lineVals[ colOrder['ALT'] ] ) )
			continue
		
		GTndx = None
		DP_AD_ndx = None
		
		lineGenoFormat = lineVals[ colOrder['FORMAT'] ].split( ':' )
		
		for ndx in range( len( lineGenoFormat ) ):
			if lineGenoFormat[ndx] == "GT":
				GTndx = ndx
			elif lineGenoFormat[ndx] == "AD" and args.filterByAD:
				DP_AD_ndx = ndx
			elif lineGenoFormat[ndx] == "DP" and not args.filterByAD:
				DP_AD_ndx = ndx
		
		if GTndx == None:
			error_msg( "GT field missing. Is this a sites file???" )
		elif DP_AD_ndx == None and not args.noFilter:
			error_msg( "You are filtering by individual sample %s field, but it wasn't found. Try filtering by %s instead?" % ( "AD" if args.filterByAD else "DP", "DP" if args.filterByAD else "AD" ) )
		
		genotypeArray = []
		
		for currSample in sampleNames:
			sampleValues = lineVals[ colOrder[currSample] ].split(":")
			
			genotype = sampleValues[GTndx][:3]
		
			if genotype == "./." or genotype == ".|.":
				genotypeArray.append( "NC" )
			else:
				if not args.noFilter:
					try:
						depthValues = sampleValues[DP_AD_ndx].split( ',' )
						depthCounts = [ int( val ) for val in depthValues ]
					except:
						genotypeArray.append( 'NC' )
						continue
				
				if args.noFilter or sum( depthCounts ) >= args.minDepth:
					keepSampleCount += 1
					genotypeArray.append( genotype_conversion( genotype ) )
				else:
					genotypeArray.append( "NC" )
		
		if keepSampleCount >= args.minCount:
			print "%s,%s,%s,%s,%s,%s" % ( lineVals[ colOrder['CHROM'] ], lineVals[ colOrder['POS'] ], lineVals[ colOrder['ID'] ], lineVals[ colOrder['REF'] ], lineVals[ colOrder['ALT'] ], ','.join( genotypeArray )  )

	VCF.close()
