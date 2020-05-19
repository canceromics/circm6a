package main;

public class InParam {
	
	private String ip_file=null;
	private String input_file=null;
	private String trim_file=null;
	private String genome_file=null;
	private String gtf_file=null;
	private String out_prefix=null;
	private String rrna_bed=null;
	private String circ_bed=null;
	private String adjust_p=null;
	private int sup_reads=2;
	private int background_size=0;
	private int window_size=25;
	private int peak_length=100;
	private int circ_length=100;
	private int circ_max_length=200000;
	private int read_dev=5;
	private double p_value=0.05;
	private boolean out_detail=false;
	private boolean pair_mode=true;
	private boolean uniq_mode=false;
	private boolean retain_test=false;
	private boolean help_flag=false;
	
	public InParam(String[] args) {
		String[] to_read={"-ip", "-input", "-trim", "-g", "-r", "-o", "-rrna", "-circ", "-adjust", "-sup", "-back", "-window", "-peak", "-cl", "-clmax", "-dev", "-pt", "-detail", "-unpair", "-uniq", "-test", "-h"};
		for (int i = 0; i < args.length; ++i) {
			if (args[i].charAt(0) == '-') {
				if (args[i].equals(to_read[0])) {
					++i;
					if (i < args.length) {
						this.ip_file = args[i];
					}
				}
				else if (args[i].equals(to_read[1])) {
					++i;
					if (i < args.length) {
						this.input_file = args[i];
					}
				}
				else if (args[i].equals(to_read[2])) {
					++i;
					if (i < args.length) {
						this.trim_file = args[i];
					}
				}
				else if (args[i].equals(to_read[3])) {
					++i;
					if (i < args.length) {
						this.genome_file = args[i];
					}
				}
				else if (args[i].equals(to_read[4])) {
					++i;
					if (i < args.length) {
						this.gtf_file = args[i];
					}
				}
				else if (args[i].equals(to_read[5])) {
					++i;
					if (i < args.length) {
						this.out_prefix = args[i];
					}
				}
				else if (args[i].equals(to_read[6])) {
					++i;
					if (i < args.length && args[i].charAt(0) != '-') {
						this.rrna_bed = args[i];
					}
					else {
						this.rrna_bed = "/rRNA.bed";
						i--;
					}
				}
				else if (args[i].equals(to_read[7])) {
					++i;
					if (i < args.length) {
						this.circ_bed = args[i];
					}
				}
				else if (args[i].equals(to_read[8])) {
					++i;
					if (i < args.length) {
						this.adjust_p = args[i];
					}
				}
				else if (args[i].equals(to_read[9])) {
					++i;
					if (i < args.length) {
						this.sup_reads = Integer.parseInt(args[i]);
						if (this.sup_reads < 0) {
							this.sup_reads = - this.sup_reads;
						}
					}
				}
				else if (args[i].equals(to_read[10])) {
					++i;
					if (i < args.length) {
						this.background_size = Integer.parseInt(args[i]);
						if (this.background_size < 0) {
							this.background_size = - this.background_size;
						}
					}
				}
				else if (args[i].equals(to_read[11])) {
					++i;
					if (i < args.length) {
						this.window_size = Integer.parseInt(args[i]);
						if (this.window_size < 0) {
							this.window_size = - this.window_size;
						}
					}
				}
				else if (args[i].equals(to_read[12])) {
					++i;
					if (i < args.length) {
						this.peak_length = Integer.parseInt(args[i]);
						if (this.peak_length < 0) {
							this.peak_length = - this.peak_length;
						}
					}
				}
				else if (args[i].equals(to_read[13])) {
					++i;
					if (i < args.length) {
						this.circ_length = Integer.parseInt(args[i]);
						if (this.circ_length < 0) {
							this.circ_length = - this.circ_length;
						}
					}
				}
				else if (args[i].equals(to_read[14])) {
					++i;
					if (i < args.length) {
						this.circ_max_length = Integer.parseInt(args[i]);
						if (this.circ_max_length < 0) {
							this.circ_max_length = - this.circ_max_length;
						}
					}
				}
				else if (args[i].equals(to_read[15])) {
					++i;
					if (i < args.length) {
						this.read_dev = Integer.parseInt(args[i]);
						if (this.read_dev < 0) {
							this.read_dev = - this.read_dev;
						}
					}
				}
				else if (args[i].equals(to_read[16])) {
					++i;
					if (i < args.length) {
						this.p_value = Double.parseDouble(args[i]);
						if (this.p_value < 0) {
							this.p_value = - this.p_value;
						}
					}
				}
				else if (args[i].equals(to_read[17])) {
					this.out_detail = true;
				}
				else if (args[i].equals(to_read[18])) {
					this.pair_mode = false;
				}
				else if (args[i].equals(to_read[19])) {
					this.uniq_mode = true;
				}
				else if (args[i].equals(to_read[20])) {
					this.retain_test = true;
				}
				else if (args[i].equals(to_read[21])) {
					this.help_flag = true;
				}
				else {
					System.out.println("Unknown parameter: " + args[i]);
				}
			}
		}
	}

	/**
	 * parameter examining
	 * @return true if pass
	 */
	public boolean completeArgs() {
		String[] for_help={"\t-ip\tIP sam/bam file searching back junctions",
				"\t-input\tINPUT sam/bam file searching back junctions, calling IP peaks with IP file)",
				"\t-trim\ttrimed sam/bam file calling peak instead of INPUT file",
				"\t-g\ta reference genome file which is the same as the mapping one",
				"\t-r\ta reference gencode file in gtf format",
				"\t-o\tprefix of out files",
				"\t-rrna\tenable rRNA remove (can specify a bed file)",
				"\t-circ\tprovide circ bed file to call peak",
				"\t-adjust\tchoose method for p-value adjusting [bon|bh](bon means Bonferroni, bh means Benjamini and Hochberg)",
				"\t-sup\tmin support number of a back junction (default is 2)",
				"\t-back\tbackground size in peak calling (0 for whole genome, default is 0)",
				"\t-window\twindow size in peak calling (default is 25)",
				"\t-peak\tmin peak length in peak calling (default is 100)",
				"\t-cl\tcutoff of distance bewteen junction points (default is 100)",
				"\t-clmax\tupper bound cutoff of distance bewteen junction points (default is 200000)",
				"\t-dev\tmax deviation permitted in searching back junctions (default is 5)",
				"\t-pt\tthreshold for p-value in peak calling (default is 0.05)",
				"\t-detail\toutput circ detail file",
				"\t-unpair\tunable pair end examination of back junction",
				"\t-uniq\tenable using only uniq mapping reads",
				"\t-h\tshow help text"
		};
		boolean out = true;
		if (this.help_flag) {
			System.out.println("For usage:");
			for (int j=0; j < for_help.length; ++j) {
				System.out.println(for_help[j]);
			}
			return false;
		}
		if (this.trim_file != null && this.input_file == null) {
			this.input_file = this.trim_file;
			this.trim_file = null;
		}
		if (this.ip_file == null && this.input_file == null) {
			this.LackParam(for_help[0]);
			System.out.println(for_help[1]);
			out = false;
		}
		if (this.out_prefix == null) {
			this.LackParam(for_help[5]);
			out = false;
		}
		if (this.window_size <= 0 || this.peak_length <= 0 || this.sup_reads <= 0 || this.background_size < 0 
				|| this.circ_length <= 0 || this.read_dev < 0 || this.p_value < 0) {
			System.out.println("ERROR: Parameter(s) exist illegal negtive number");
			out = false;
		}
		return out;
	}
	
	public void LackParam(String s) {
		System.out.println("ERROR: Missing parameter:" + s);
	}
	
	public String getIp_file() {
		return ip_file;
	}

	public String getInput_file() {
		return input_file;
	}

	public String getGenome_file() {
		return genome_file;
	}

	public String getGtf_file() {
		return gtf_file;
	}

	public String getOut_prefix() {
		return out_prefix;
	}

	public String getRrna_bed() {
		return rrna_bed;
	}
	
	public String getCirc_bed() {
		return circ_bed;
	}

	public int getSup_reads() {
		return sup_reads;
	}

	public int getBackground_size() {
		return background_size;
	}

	public int getWindow_size() {
		return window_size;
	}

	public int getPeak_length() {
		return peak_length;
	}

	public int getRead_dev() {
		return read_dev;
	}

	public double getP_value() {
		return p_value;
	}

	public boolean isOut_detail() {
		return out_detail;
	}

	public boolean isPair_mode() {
		return pair_mode;
	}

	public boolean isUniq_mode() {
		return uniq_mode;
	}
	
	public boolean isRetain_test() {
		return retain_test;
	}
	
	public String getTrim_file() {
		return trim_file;
	}

	public int getCirc_length() {
		return circ_length;
	}
	
	public int getCirc_max_length() {
		return circ_max_length;
	}

	public String getAdjsut_p() {
		return adjust_p;
	}
}