#!/bin/env/ruby

#----------------------------------------------------------------------------------------
# uce2speciestree
UCE2SPECIESTREEVER = "0.4.0"
# Michael G. Campana, 2019-2025  
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require 'optparse'
require 'ostruct'
$filecount = 0 # Max number of files
$samples = {} # Hash of samples
$loci = {} # Hash of loci and lengths
#----------------------------------------------------------------------------------------
class Sample
	attr_accessor :name, :loci, :concat
	def initialize(name, loci = {}, concat = "")
		@name = name # Sequence name
		@loci = loci # Hash of loci and sequences
		@concat = concat # Final concatenated sequence
	end
	
end
#----------------------------------------------------------------------------------------
def get_files
	@file_names = "OriginalLocus\tRenamedLocus\n"
	Dir.foreach(File.expand_path($options.indir)) do |file|
		unless File.directory?(file)
			if file[-4..-1] == ".nex" or file[-6..-1] == ".nexus"
				$filecount += 1
				$loci[file] = ""
				@file_names << $options.indir + "/" + file + "\t" + $options.outdir + "/locus." + $filecount.to_s + ".fa\n"
				convert_nexus(file)
			end
		end
	end
	File.open($options.outdir + "/locusnames.tsv", 'w') do |write|
		write.puts @file_names
	end
end
#----------------------------------------------------------------------------------------
def convert_nexus(file)
	@out = ""
	start = false
	File.open(File.expand_path($options.indir) + "/" + file, 'r') do |getseq|
		while line = getseq.gets
			start = false if line == ";\n"
			if !start && line.include?("nchar=")
				line_arr = line.split("nchar=")
				$loci[file] = line_arr[1][0..-2].to_i
			elsif start
				unless  line == "\n"
					line_arr = line.split(" ")
					make_loci(file,line_arr[0],convert_seq(line_arr[-1]))
				end
			end
			start = true if line == "matrix\n"
		end
	end
	for sample in $samples.keys
		unless $samples[sample].loci[file].nil?
			@out << ">" + sample + "\n" #Append header to outline
			@out <<  $samples[sample].loci[file] + "\n" # Append sequence to outline
		end
	end
	File.open(File.expand_path($options.outdir) + "/locus." + $filecount.to_s + ".fa", 'w') do |write|
		write.puts @out
	end
end
#----------------------------------------------------------------------------------------
def make_loci(locus,name,sequence)
	if !$samples.keys.include?(name)
		$samples[name] = Sample.new(name, { locus => sequence })
	elsif !$samples[name].loci.keys.include?(locus)
		$samples[name].loci[locus] = sequence
	else
		$samples[name].loci[locus] << sequence
	end
end
#----------------------------------------------------------------------------------------
def make_paup
	@total_data_set = 0
	@charset = "begin sets;\n"
	@charpar = "\n    charpartition combined = "
	for locus in $loci.keys
		@total_data_set += $loci[locus]
		@charset << "    charset " + locus.delete("-") + " = " + (@total_data_set - $loci[locus] + 1).to_s + "-" + @total_data_set.to_s + ";\n"
		@charpar << locus.delete("-") + ":" + (@total_data_set - $loci[locus] + 1).to_s + "-" + @total_data_set.to_s + ", "
		for sample in $samples.keys # Find missing data
			if $samples[sample].loci[locus].nil?
				$samples[sample].loci[locus] = "?" * $loci[locus]
			end
			$samples[sample].concat << $samples[sample].loci[locus]
		end
	end
	@charpar[-2..-1] = ";\nend;\n" # Replace final locus
	@header = "#NEXUS\nbegin data;\n    dimensions ntax=#{$samples.keys.size} nchar=#{@total_data_set};\n    format datatype=nucleotide missing=? gap=-;\nmatrix\n"
	for sample in $samples.keys
		@header << "    " + sample + "    " + $samples[sample].concat + "\n"
	end
	@header << ";\nend;\n\n" + @charset + @charpar + "\n\nbegin paup;\n    svdq evalq=all bootstrap;"
	@header << "\n    Outgroup #{$options.paupout};" unless $options.paupout.nil?
	@header << "\n    SaveTrees file=#{$options.paupsave} format=#{$options.paupformat};\n    quit;\nend;\n"
	File.open($options.outdir + "/uce2svdquartets.nex", 'w') do |write|
		write.puts @header
	end
end
#----------------------------------------------------------------------------------------
def convert_seq(fasta)
	for i in 0 ... fasta.length
		break if fasta[i].chr != "-"
		fasta[i] = "?"
	end
	if fasta.length > 1
		i = fasta.length - 1 # Get end of sequence
		while fasta[i].chr == "-"
			fasta[i] = "?"
			i -= 1
		end
	end
	return fasta
end
#----------------------------------------------------------------------------------------
def write_array_qsub
	header = "#!/bin/sh\n#$ -S /bin/sh\n#$ -q #{$options.queue}\n#$ -l mres=#{$options.mres},h_data=#{$options.h_data},h_vmem=#{$options.vmem}"
	header << ",himem" if $options.himem
	header << ",lopri" if $options.lopri
	header << "\n#$ -j y\n#$ -cwd\n#$ -N uce2speciestree\n"
	header << "#$ -m bea\n#$ -M #{$options.email}\n" if $options.email != ""
	header << "module load bio/raxml\nraxmlHPC-SSE3 -s #{$options.outdir}/locus.$SGE_TASK_ID.fa -w #{File.expand_path($options.outdir)} -n locus.$SGE_TASK_ID.tree -m #{$options.submodel} -p #{$options.raxp} -N #{$options.raxn}"
	if $options.ml
		header << " -f a -x #{$options.raxb[0]}"
	elsif $options.raxb[1]
		header << " -b #{$options.raxb[0]}"
	end
	header <<  " #{$options.raxml}\n"
	File.open($options.outdir + '/uce2speciestree.job', 'w') do |qsub|
		qsub.puts header
	end
end
#----------------------------------------------------------------------------------------
def execute_qsub_array
	system("qsub -t 1-#{$filecount}:1  #{$options.outdir + '/uce2speciestree.job'}")
end
#----------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.indir = "./" # Input directory
		args.outdir = "./" # Output directory
		args.mres = "1G" # Memory reserved
		args.h_data = "1G" # h_data per CPU
		args.vmem = "1G" # Virtual memory
		args.himem = false # High memory flag
		args.lopri = false # Low priority falg
		args.queue = "" # Qsub queue
		args.email = "" # Email address to notify
		args.submodel = "GTRGAMMA" # Substitution model to pass to RAxML
		args.raxp = srand # Random seed for RAxML
		args.raxb = [srand,false] # Random seed for RAxML bootstrap
		args.ml = false # Perform ML search in RAxML
		args.raxn = 100 # Number of RAxML bootstrap replicates
		args.raxml = "" # String of parameters to pass to RAxML
		args.paupsave = "uce2svdquartets.tre" # PAUP* saved tree filename
		args.paupformat = "Newick" # PAUP* saved tree file format
		args.paupout = nil # PAUP* outgroup
		opt_parser = OptionParser.new do |opts|
			opts.banner = "Command-line usage: ruby uce2speciestree.rb [options]"
			opts.separator ""
			opts.separator "Input/Output Options:"
			opts.on("-i","--input [FILE]", String, "Input Directory (Default is current directory)") do |indir|
				args.indir = File.expand_path(indir) if indir != nil
			end
			opts.on("-o","--output [FILE]", String, "Output Directory (Default is current directory)") do |outdir|
				args.outdir = File.expand_path(outdir)if outdir != nil
			end
			opts.separator ""
			opts.separator "RAxML options:"
			opts.on("-d", "--submodel [VALUE]", String, "RAxML substitution model (Default = GTRGAMMA)") do |submodel|
				args.submodel = submodel if submodel != nil
			end
			opts.on("-p", "--raxp [VALUE]", Integer, "RAxML random number seed (Default uses system entropy)") do |raxp|
				args.raxp = raxp if raxp != nil
			end
			opts.on("-b", "--raxb [VALUE]", Integer, "Bootstrap RAxML and set random number seed. Omitting value for flag uses system entropy.") do |raxb|
				args.raxb[1] = true
				args.raxb[0] = raxb if raxb != nil
			end
			opts.on("-M", "--ML", "Perform ML search after bootstrapping (converts -b to -x and adds -f a to RAxML)") do |ml|
				args.ml = true
				args.raxb[1] = true
			end
			opts.on("-N", "--raxn [VALUE]", Integer, "Number of RAxML alternative runs (Default = 100)") do |raxn|
				args.raxn = raxn if raxn != nil
			end
			opts.on("-r[VALUE]", String, "String of RAxML optional parameters. Enclose string in '' and leave no space between the string and the flag.") do |raxml|
				args.raxml = raxml if raxml != nil
			end
			opts.separator ""
			opts.separator "Hydra options:"
			opts.on("-q", "--queue [VALUE]", String, "Qsub queue to use") do |queue|
				args.queue = queue if queue != nil
			end
			opts.on("-m", "--mres [VALUE]", String, "Total reserved memory (Default = 1G)") do |mres|
				args.mres = mres if mres != nil
			end
			opts.on("-t", "--h_data [VALUE]", String, "h_data per CPU (Default = 1G)") do |h_data|
				args.h_data = h_data if h_data != nil
			end
			opts.on("-V", "--vmem [VALUE]", String, "Total reserved virtual memory (Default = 1G)") do |vmem|
				args.vmem = vmem if vmem != nil
			end
			opts.on("-H", "--himem", "Use high-memory queue (Default is false)") do |himem|
				args.himem = true
			end
			opts.on("-L", "--lopri", "Use low priority queue (Default is false)") do |lopri|
				args.lopri = true
			end
			opts.on("-e", "--email [VALUE]", String, "E-mail address to notify") do |email|
				args.email = email if email != nil
			end
			opts.separator ""
			opts.separator "PAUP* options:"
			opts.on("-s", "--paupsave [VALUE]", String, "Saved PAUP* tree file name (Default = uce2svdquartets.tre)") do |paupsave|
				args.paupsave = paupsave if paupsave != nil
			end
			opts.on("-F", "--paupformat [VALUE]", String, "Saved PAUP* tree file format (Default = Newick)") do |paupformat|
				args.paupformat = paupformat if paupformat != nil
			end
			opts.on("-O", "--paupout [VALUE]", String, "PAUP* outgroup (Default = unrooted)") do |paupout|
				args.paupout = paupout if paupout != nil
			end
			opts.separator ""
			opts.separator "General information:"
			opts.on("-v", "--version", "Show pipeline program version(s)") do
				puts "uce2speciestree.rb " + UCE2SPECIESTREEVER + "\n"
				exit
			end
			opts.on_tail("-h","--help", "Show help") do
				puts opts
				exit
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#--------------------------------------------------------
ARGV[0] ||= "-h"
$options = Parser.parse(ARGV)
Dir.mkdir($options.outdir) if !FileTest.directory?($options.outdir)
get_files
make_paup
write_array_qsub
execute_qsub_array