# huntcircRNA
this is a practice till now
For search BSJ in Sam/Bam file(s)
To use:
java -jar -ip IP_file <-input input_file> <-trim trim_file> -o out_prefix <options>
		-ip		[String]	IP sam/bam file searching back junctions,
		-input	[String]	INPUT sam/bam file searching back junctions, calling IP peaks with IP file)
		-trim	[String]	trimed sam/bam file calling peak instead of INPUT file,
		-o		[String]	prefix of out files,
options:-g		[String]	a reference genome file which is the same as the mapping one,
		-r		[String]	a reference gencode file in gtf format,
		-rrna	<String>	enable rRNA remove (can specify a bed file),
		-circ	[String]	provide circ bed file to call peak,
		-sup	[int]		min support number of a back junction (default is 1),
		-back	[int]		background size in peak calling (0 for whole genome, default is 0),
		-window	[int]		window size in peak calling (default is 25),
		-peak	[int]		min peak length in peak calling (default is 100),
		-cl		[int]		cutoff of distance bewteen junction points (default is 100),
		-dev	[int]		max deviation permitted in searching back junctions (default is 5),
		-pt		[float]		threshold for p-value in peak calling (default is 0.05),
		-detail				output circ detail file,
		-pair				enable pair end examination of back junction,
		-h					show help text