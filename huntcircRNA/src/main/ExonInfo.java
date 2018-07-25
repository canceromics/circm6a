package main;

public class ExonInfo {
	private String base_seq=null;
	private String chr_symbol=null;
	private int chr_num=0;
	private int start_position=0;
	private int end_position=0;
	private int read_count=0;
	private ExonInfo next=null;
	
	String getBase_seq() {
		return base_seq;
	}

	void setBase_seq(String base_seq) {
		this.base_seq = base_seq;
	}

	String getChr_symbol() {
		return chr_symbol;
	}
	
	void setChr_symbol(String chr_symbol) {
		this.chr_symbol = chr_symbol;
		this.chr_num = chrSymbolToNum(chr_symbol);
	}
	
	int getChr_num() {
		return chr_num;
	}
	
	void setChr_num(int chr_num) {
		this.chr_num = chr_num;
		this.chr_symbol = chrNumToSymbol(chr_num);
	}
	
	int getStart_position() {
		return start_position;
	}
	
	void setStart_position(int start_position) {
		this.start_position = start_position;
	}
	
	int getEnd_position() {
		return end_position;
	}
	
	void setEnd_position(int end_position) {
		this.end_position = end_position;
	}
	
	int getRead_count() {
		return read_count;
	}

	void increaseRead_count() {
		this.read_count++;
	}
	
	void setRead_count(int read_count) {
		this.read_count = read_count;
	}

	public ExonInfo getNext() {
		return next;
	}

	public void setNext(ExonInfo next) {
		this.next = next;
	}

	/**
	 * transform "chr1~22XYM" into 0~24
	 * @param chr_symbol chromosome symbol like "chr1"; 
	 * @return index of the chromosome
	 */
	public static int chrSymbolToNum(String chr_symbol) {
		int out = 0;
		if ((chr_symbol.length() == 4 || chr_symbol.length() == 5) && chr_symbol.matches("chr[0-9XYM][0-9]{0,1}")) {
			String temp = chr_symbol.substring(3);
			if (temp.equals("X")) {
				out = 23;
			}
			else if (temp.equals("Y")) {
				out = 24;
			}
			else if (temp.equals("M")) {
				out = 25;
			}
			else {
				out = Integer.parseInt(temp);
			}
		}
		--out;
		return out;
	}
	
	/**
	 * transform 0~24 into "chr1~22XYM"
	 * @param chr_num index of the chromosome;
	 * @return chromosome symbol
	 */
	public static String chrNumToSymbol(int chr_num) {
		String out = null;
		++chr_num;
		if (chr_num > 25 || chr_num <= 0) {
			 System.out.println("Undefine CHR");
		}
		else if (chr_num==23) {
			out = "chrX";
		}
		else if (chr_num==24) {
			out = "chrY";
		}
		else if (chr_num==25) {
			out = "chrM";
		}
		else {
			out = "chr" + chr_num;
		}
		return out;
	}
	
	/**
	 * calculate end position according to start position and cigar string
	 * @param start start position
	 * @param cigar cigar string
	 * @return end position
	 */
	public static int getEnd(int start, String cigar) {
		int out = start - 1;
		for (int i = 0; i < cigar.length(); ++i) {
			char c = cigar.charAt(i);
			int seg = 0;
			while (c <= '9' && c >= '0') {
				seg = seg * 10 + c - 48;// '0' == 48
				++i;
				c = cigar.charAt(i);
			}
			if (c == 'M' || c == 'N' || c == 'D') {
				out += seg;
			}
			seg = 0;
		}
		return out;
	}
	
}
