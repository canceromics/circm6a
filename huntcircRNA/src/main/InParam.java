package main;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

public class InParam {
	
	private static InParam param = new InParam();
	private String ip_file=null;
	private String input_file=null;
	private String trim_file=null;
	private String genome_file=null;
	private String gtf_file=null;
	private String out_prefix=null;
	private String rrna_bed=null;
	private String circ_bed=null;
	private String adjust_p=null;
	private String combine_p="log";
	private int sup_reads=2;
	private int ip_sup=2;
	private int window_size=25;
	private int peak_length=100;
	private int circ_length=100;
	private int circ_max_length=200000;
	private int read_dev=8;
	private int map_quality=10;
	private double p_value=0.05;
	private double propor=1.0;
	private boolean peak_BSJ=false;
	private boolean out_detail=false;
	private boolean peak_quali=false;
	private boolean XA_mode=false;
	private boolean uniq_mode=false;
//	private boolean retain_test=false;
	private boolean help_flag=false;
	
	
	private InParam() {
	}
	
	public static InParam getParams() {
		return param;
	}
	
	public void SetParams(String[] args) {
		for (int i = 0; i < args.length; ++i) {
			if (args[i].charAt(0) == '-') {
				switch (args[i].toLowerCase()) {
				case "-ip":
					param.ip_file = ++i < args.length ? args[i] : param.ip_file;
					break;

				case "-input":
					param.input_file = ++i < args.length ? args[i] : param.input_file;
					break;
					
				case "-trim":
					param.trim_file = ++i < args.length ? args[i] : param.trim_file;
					break;
					
				case "-g":
					param.genome_file = ++i < args.length ? args[i] : param.genome_file;
					break;
					
				case "-r":
					param.gtf_file = ++i < args.length ? args[i] : param.gtf_file;
					break;
					
				case "-o":
					param.out_prefix = ++i < args.length ? args[i] : param.out_prefix;
					break;
					
				case "-rrna":
					param.rrna_bed = i + 1 < args.length && args[i + 1].charAt(0) != '-' ? args[++i] : "/rRNA.bed";
					break;
					
				case "-circ":
					param.circ_bed = ++i < args.length ? args[i] : param.circ_bed;
					break;
					
				case "-adjust":
					param.adjust_p = ++i < args.length ? args[i] : param.adjust_p;
					break;
					
				case "-combine":
					param.combine_p = ++i < args.length ? args[i] : param.adjust_p;
					break;
					
				case "-sup":
					param.sup_reads = ++i < args.length ? Math.abs(Integer.parseInt(args[i])) : param.sup_reads;
					break;
					
				case "-ipsup":
					param.ip_sup = ++i < args.length ? Math.abs(Integer.parseInt(args[i])) : param.ip_sup;
					break;
					
				case "-ippro":
					param.propor = ++i < args.length ? Math.abs(Double.parseDouble(args[i])) : param.propor;
					break;
					
				case "-window":
					param.window_size = ++i < args.length ? Math.abs(Integer.parseInt(args[i])) : param.window_size;
					break;
					
				case "-peak":
					param.peak_length = ++i < args.length ? Math.abs(Integer.parseInt(args[i])) : param.peak_length;
					break;
					
				case "-cl":
					param.circ_length = ++i < args.length ? Math.abs(Integer.parseInt(args[i])) : param.circ_length;
					break;
					
				case "-clmax":
					param.circ_max_length = ++i < args.length ? Math.abs(Integer.parseInt(args[i])) : param.circ_max_length;
					break;
					
				case "-mapq":
					param.map_quality = ++i < args.length ? Math.abs(Integer.parseInt(args[i])) : param.map_quality;
					break;
					
				case "-dev":
					param.read_dev = ++i < args.length ? Math.abs(Integer.parseInt(args[i])) : param.read_dev;
					break;
					
				case "-pt":
					param.p_value = ++i < args.length ? Math.abs(Double.parseDouble(args[i])) : param.p_value;
					break;
					
				case "-detail":
					param.out_detail = true;
					break;
					
				case "-xa":
					param.XA_mode = true;
					break;
					
				case "-uniq":
					param.uniq_mode = true;
					break;
					
				case "-peakq":
					param.peak_quali = true;
					break;
					
				case "-peakb":
					param.peak_BSJ = true;
					break;
//					
//				case "-test":
//					param.retain_test = true;
//					break;
					
				case "-h":
					param.help_flag = true;
					break;
					
				default:
					System.out.println("Unknown parameter: " + args[i]);
					break;
				}
			}
			else {
				System.out.println("Unknown parameter: " + args[i]);
			}
		}
	}

	/**
	 * parameter examining
	 * @return true if pass
	 */
	public boolean completeArgs() {
		boolean out = true;
		if (this.help_flag) {
			try (BufferedReader reader = new BufferedReader(new InputStreamReader(RemoverRNA.class.getResourceAsStream("/HelpDocument"), "UTF-8"))){
				String line = null;
				while ((line = reader.readLine()) != null) {
					String[] cols = line.split("\t");
					if (cols.length == 3) {
						System.out.printf("\t%-24s%s\n", cols[1], cols[2]);
					}
					else {
						System.out.println(line);
					}
					
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
			return false;
		}
		if (this.trim_file != null && this.input_file == null) {
			this.input_file = this.trim_file;
			this.trim_file = null;
		}
		if (this.input_file == null) {
			this.LackParam("<input.bam>");
			out = false;
		}
		if (this.out_prefix == null) {
			this.LackParam("<path/out_prefix>");
			out = false;
		}
		if (this.window_size <= 0 || this.peak_length <= 0 || this.sup_reads <= 0 || this.ip_sup < 0
				|| this.circ_length <= 0 || this.read_dev < 0 || this.p_value < 0) {
			System.out.println("ERROR: Parameter(s) exist illegal negtive number or 0");
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

	public int getIp_sup() {
		return ip_sup;
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
	
	public boolean isUniq_mode() {
		return uniq_mode;
	}
	
	public String getCombine_p() {
		return combine_p;
	}
//	
//	public boolean isRetain_test() {
//		return retain_test;
//	}
	
	public String getTrim_file() {
		return trim_file;
	}

	public int getCirc_length() {
		return circ_length;
	}
	
	public int getCirc_max_length() {
		return circ_max_length;
	}

	public int getMapQuality() {
		return map_quality;
	}

	public double getPropor() {
		return propor;
	}
	
	public String getAdjsut_p() {
		return adjust_p;
	}
	
	public boolean isXA_mode() {
		return XA_mode;
	}
	
	public boolean isPeakQuality() {
		return peak_quali;
	}
	
	public boolean isPeakBSJ() {
		return peak_BSJ;
	}
	
	void setPropor(double propor) {
		this.propor = propor;
	}
	
	
}