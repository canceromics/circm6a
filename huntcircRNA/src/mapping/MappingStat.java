package mapping;

import main.InParam;

public class MappingStat {
	
	private static long ip_reads = 0;
	private static long mapping_ip_reads = 0;
	private static long input_reads = 0;
	private static long mapping_input_reads = 0;
	
	private int reads = 0;
	private int map_reads = 0;
	private int reads_ip = 0;
	private int map_reads_ip = 0;
	private int circ = 0;
	private int circ_ip = 0;
	
	int incReads(boolean ip_flag, int mapq) {
		return mapq >= InParam.getParams().getMapQuality() ? (ip_flag ? ++map_reads_ip + ++reads_ip : ++map_reads + ++reads)
				: (ip_flag ? ++reads_ip : ++reads);
	}
	
	void resetReads(boolean ip_flag) {
		if (ip_flag) {
			reads_ip = 0;
			map_reads_ip = 0;
			circ_ip = 0;
		}
		else {
			reads = 0;
			map_reads = 0;
			circ = 0;
		}
	}
	
	int incCircReads(boolean ip_flag) {
		return ip_flag ? ++circ_ip : ++circ;
	}
	
	public int getReads(boolean ip) {
		return ip ? reads_ip : reads;
	}
	
	public int getReads() {
		return reads;
	}
	
	public int getMapReads() {
		return map_reads;
	}
	
	public int getReadsIP() {
		return reads_ip;
	}
	
	public int getMapReadsIP() {
		return map_reads_ip;
	}
	
	public int getCircReads() {
		return circ;
	}
	
	public int getCircReadsIP() {
		return circ_ip;
	}
	
	public static long getIpReads() {
		return ip_reads;
	}
	
	public static long getMappingIpReads() {
		return mapping_ip_reads;
	}
	
	public static long getInputReads() {
		return input_reads;
	}
	
	public static long getMappingInputReads() {
		return mapping_input_reads;
	}
	
	static long incReadNum(boolean ip_flag, int mapq) {
		return mapq >= InParam.getParams().getMapQuality() ? 
			(ip_flag ? ++mapping_ip_reads + ++ip_reads : ++mapping_input_reads + ++input_reads) : (ip_flag ? ++ip_reads : ++input_reads);
	}
	
	static void resetReadNum(boolean ip_flag) {
		if (ip_flag) {
			mapping_ip_reads = 0;
			ip_reads = 0;
		}
		else {
			mapping_input_reads = 0;
			input_reads = 0;
		}
	}
	
}
