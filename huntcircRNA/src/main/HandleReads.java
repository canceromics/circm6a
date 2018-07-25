package main;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map.Entry;

import htsjdk.samtools.util.IntervalTree;


public class HandleReads {
	
	private int read_lenth = 0;
	private int seg_dev = 8;
	private class Segment{
		private int chr_num=0;
		private int pos_start=0;
		private int pos_end=0;
		private boolean positive;
	}
	
	/**
	 * find pcc signal in alignments with the same id
	 * @param reads first alignments with this id;
	 * @param length count first alignments;
	 * @param out junctions found(output);
	 * @param reads2 last alignments with this id;
	 * @param length2 count first alignments;
	 * @param ip_flag whether these alignments in IP file;
	 * @param itree tree to store alignments;
	 * @param circ_length min length between two junction site;
	 * @param r removing rRNA;
	 * @param pair_flag whether enable pair examine;
	 * @return removed rRNA alignments;
	 */
	public int filtPCC(String[] reads, int length, ArrayList<HashMap<String, JuncInfo>> out,String[] reads2, int length2, boolean ip_flag, ArrayList<IntervalTree<ReadInfo>> itree, int circ_length, RemoverRNA r, boolean pair_flag, boolean add_flag, boolean uniq_map) {
		int rrna = 0;
		if (length >= 1) {
			HashMap<Integer, ArrayList<Segment>> seg_map = new HashMap<>();
			HashMap<Integer, ArrayList<Segment>> seg_map2 = new HashMap<>();
			ArrayList<Segment> examin_list1 = new ArrayList<>();
			ArrayList<Segment> examin_list2 = new ArrayList<>();
			String id = null;
			for (int i = 0; i < length; ++i) {
				String[] cols = reads[i].split("\t");
				id = cols[0];
				int chr_num = ExonInfo.chrSymbolToNum(cols[2]);
				int xa_index = reads[i].indexOf("XA:Z:");
				if (chr_num >= 0 && chr_num < out.size() && (!uniq_map || xa_index == -1)) {
					this.addSegment(seg_map, chr_num, cols[5], Integer.parseInt(cols[3]), (Integer.parseInt(cols[1]) >>> 4 & 1) == 0, r, examin_list1);
				}
				if (!uniq_map && xa_index != -1) {
					String xa = reads[i].substring(xa_index);
					int last_col = xa.indexOf("\t");
					if (last_col == -1) {
						xa = xa.substring(xa.lastIndexOf(':') + 1);
					}
					else {
						xa = xa.substring(0, last_col);
						xa = xa.substring(xa.lastIndexOf(':') + 1);
					}
					this.putXAstate(seg_map, xa, r);
				}
			}
			for (int i = 0; i < length2; ++i) {
				String[] cols = reads2[i].split("\t");
				id = cols[0];
				int chr_num = ExonInfo.chrSymbolToNum(cols[2]);
				int xa_index = reads2[i].indexOf("XA:Z:");
				if (chr_num >= 0 && chr_num < out.size() && (!uniq_map || xa_index == -1)) {
					this.addSegment(seg_map, chr_num, cols[5], Integer.parseInt(cols[3]), (Integer.parseInt(cols[1]) >>> 4 & 1) == 0, r, examin_list2);
				}
				if (!uniq_map && xa_index != -1) {
					String xa = reads2[i].substring(xa_index);
					int last_col = xa.indexOf("\t");
					if (last_col == -1) {
						xa = xa.substring(xa.lastIndexOf(':') + 1);
					}
					else {
						xa = xa.substring(0, last_col);
						xa = xa.substring(xa.lastIndexOf(':') + 1);
					}
					this.putXAstate(seg_map, xa, r);
				}
			}
			if (seg_map.size() >= 1) {
				ArrayList<Segment> clip_list = new ArrayList<>();
				for (Integer key : seg_map.keySet()) {
					this.getAllClip(id, key, this.seg_dev, seg_map, clip_list, itree, examin_list2);
				}
				for (int i = 0; i < clip_list.size(); ++i) {
					Segment clip_seg = clip_list.get(i);
					if (clip_seg.pos_end - clip_seg.pos_start > circ_length) {
						String key = clip_seg.pos_start + "\t" + clip_seg.pos_end;
						HashMap <String, JuncInfo> the_chr = out.get(clip_seg.chr_num);
						if (!the_chr.containsKey(key)) {
							JuncInfo the_junc = null;
							the_junc = new JuncInfo(clip_seg.pos_start, clip_seg.pos_end);
							if (add_flag && (!pair_flag || examin_list2.size() == 0 || this.pairEndPCC(the_junc, examin_list2, clip_seg.chr_num))) {
								the_junc.addReadID(id, ip_flag);
								the_junc.setJunc_flag(ip_flag);
								the_chr.put(key, the_junc);
							}
						}
						else {
							JuncInfo the_junc = the_chr.get(key);
							if (!pair_flag || examin_list2.size() == 0 || this.pairEndPCC(the_junc, examin_list2, clip_seg.chr_num)) {
								the_junc.addReadID(id, ip_flag);
								the_junc.setJunc_flag(ip_flag);
							}
						}
					}
				}
			}
			if (seg_map2.size() >= 1) {
				ArrayList<Segment> clip_list = new ArrayList<>();
				for (Integer key : seg_map2.keySet()) {
					this.getAllClip(id, key, this.seg_dev, seg_map2, clip_list, itree, examin_list1);
				}
				for (int i = 0; i < clip_list.size(); ++i) {
					Segment clip_seg = clip_list.get(i);
					if (clip_seg.pos_end - clip_seg.pos_start > circ_length) {
						String key = clip_seg.pos_start + "\t" + clip_seg.pos_end;
						HashMap <String, JuncInfo> the_chr = out.get(clip_seg.chr_num);
						if (!the_chr.containsKey(key)) {
							JuncInfo the_junc = null;
							the_junc = new JuncInfo(clip_seg.pos_start, clip_seg.pos_end);
							if (add_flag && (!pair_flag || examin_list1.size() == 0 || this.pairEndPCC(the_junc, examin_list1, clip_seg.chr_num))) {
								the_junc.addReadID(id, ip_flag);
								the_junc.setJunc_flag(ip_flag);
								the_chr.put(key, the_junc);
							}
	//						reads.add(0, this.getchrSym(clip_seg.chr_num) + "\t" + key);
	//						FileRW.fileAppend("init.info", reads);
						}
						else {
							JuncInfo the_junc = the_chr.get(key);
							if (!pair_flag || examin_list1.size() == 0 || this.pairEndPCC(the_junc, examin_list1, clip_seg.chr_num)) {
								the_junc.addReadID(id, ip_flag);
								the_junc.setJunc_flag(ip_flag);
							}
						}
					}
				}
			}
		}
		return rrna;
	}
	
	/**
	 * decide whether the junction pair end alignment is bewteen junction sites
	 * @param junc the junction;
	 * @param segs pair-end alignments;
	 * @param chr which chr find the junction; 
	 * @return true if pass the examine
	 */
	public boolean pairEndPCC(JuncInfo junc, ArrayList<Segment> segs, int chr) {
		boolean out = false;
		for (int i = 0; i < segs.size(); ++i) {
			Segment seg = segs.get(i);
			if (seg.chr_num == chr && seg.pos_start >= junc.getSP() - seg_dev && seg.pos_end <= junc.getEP() + seg_dev) {
				out |= true;
				break;
			}
		}
		return out;
	}
	
	/**
	 * put read into segments
	 * @param seg_map all segs mapping;
	 * @param chr_num chromosome index;
	 * @param cigar_string string of cigar;
	 * @param pos start position;
	 * @param positive whether the alignment is positive mapping;
	 * @param r removing rRNA;
	 * @param pair_examin add the segment to examine(output);
	 * @return key of this segment
	 */
	public int addSegment(HashMap<Integer, ArrayList<Segment>> seg_map, int chr_num, String cigar_string, int pos, boolean positive, RemoverRNA r, ArrayList<Segment> pair_examin) {
		ArrayList<Segment> seg_list = null;
		int[] cigar = getSegment(cigar_string, positive);
		int seg_key = cigar[2] * this.read_lenth + cigar[0];
		if (seg_map.containsKey(seg_key)) {
			seg_list = seg_map.get(seg_key);
		}
		else {
			seg_list = new ArrayList<>();
		}
		Segment seg = new Segment();
		seg.chr_num = chr_num;
		seg.positive = positive;
		seg.pos_start = pos;
		seg.pos_end = pos + cigar[1];
		pair_examin.add(seg);
		if (!r.isrRNA(chr_num, seg.pos_start, seg.pos_end)) {
			seg_list.add(seg);
			seg_map.put(seg_key, seg_list);
		}
		return seg_key;
	}
	
	/**
	 * calculate cigar information
	 * @param cigar cigar characters
	 * @param positive whether + strand aligment
	 * @return list of 3 numbers:
	 * 			[0] index of the aligment end in this read
	 * 			[1] aligment numbers in reference chromosome
	 * 			[2] index of the aligment start in this read
	 */
	public int[] getSegment(String cigar, boolean positive) {
		int[] out = {0, 0, 0};
		this.read_lenth = 0;
		int start = -1;
		for (int i = 0; i < cigar.length(); ++i) {
			char c = cigar.charAt(i);
			int seg = 0;
			while (c <= '9' && c >= '0') {
				seg = seg * 10 + c - 48;// '0' == 48
				++i;
				c = cigar.charAt(i);
			}
			if (c == 'M') {
				if (start == -1) {
					start = out[0];
				}
				out[1] += seg;
				out[0] += seg;
				this.read_lenth += seg;
			}
			else if (c == 'I') {
				out[0] += seg;
				this.read_lenth += seg;
			}
			else if (c == 'N' || c == 'D') {
				out[1] += seg;
			}
			else if (c == 'S' || c == 'H') {
				if (start == -1) {
					out[0] += seg;
				}
				this.read_lenth += seg;
			}
			seg = 0;
		}
		out[0]--;
		out[2] = start;
		if (! positive) {
			out[0] = this.read_lenth - out[0] - 1;
			out[2] = this.read_lenth - out[2] - 1;
			out[0] ^= out[2];
			out[2] ^= out[0];
			out[0] ^= out[2];
		}
		out[1]--;
		return out;
	}
	
	/**
	 * fix offset of postion to the cutters nearest to mapping offset
	 * @param Dev deviation permitted to find cutters
	 * @param mapOff offset of the mapping
	 * @param cutter cutter strings
	 * @return the fixed offset
	 */
	public int fixPosOff(String Dev, int mapOff, String cutter){
		int out = -Dev.length();
		int off = -1;
		while ((off = Dev.indexOf(cutter, off + 1)) != -1){
			if (off == mapOff){
				out = off;
				break;
			}
			else {
				if (Math.abs(out - mapOff) > Math.abs(off - mapOff)) {
					out = off;
				}
			}
		}
		return out;
	}
	
	/**
	 * fix offset of postion to the cutters nearest to mapping offset
	 * @param AGDev string of deviation permitted to find AG cutters
	 * @param GTDev string of deviation permitted to find GT cutters
	 * @param AGindex index of former mapping
	 * @param GTindex index of former mapping
	 * @param AGcutter AG cutter strings
	 * @param GTcutter GT cutter strings
	 * @return vector of 2 numbers:
	 * 			[0] AGindex fixed, if Integer.MIN_VALUE appeared means cannot fix
	 * 			[1] GTindex fixed, if Integer.MIN_VALUE appeared means cannot fix
	 */
	public int[] fixPosOff(String AGDev, String GTDev, int AGindex, int GTindex, String[] AGcutter, String[] GTcutter){
		int[] out = {Integer.MIN_VALUE, Integer.MIN_VALUE};
		for (int i = 0; i < AGcutter.length; ++i){
			int off = -1;
			if (AGDev != null && GTDev != null) {
				off = this.fixPosOff(AGDev, AGindex - 3, AGcutter[i]);
				if (off > -AGDev.length()) {
					out[0] = off;
					out[0] = 3 + out[0] - AGindex;
				}
				off = this.fixPosOff(GTDev, GTindex, GTcutter[i]);
				if (off > -GTDev.length()) {
					out[1] = off;
					out[1] = out[1] - GTindex;
				}
			}
			else if (AGDev != null) {
				off = this.fixPosOff(AGDev, AGindex - 3, AGcutter[i]);
				if (off > -AGDev.length()) {
					out[0] = off;
					out[0] = 3 + out[0] - AGindex;
				}
				out[1] = 0;
			}
			else if (GTDev != null) {
				off = this.fixPosOff(GTDev, GTindex, GTcutter[i]);
				if (off > -GTDev.length()) {
					out[1] = off;
					out[1] = out[1] - GTindex;
				}
				out[0] = 0;
			}
			else {
				out[0] = 0;
				out[1] = 0;
			}
			if (out[0] == Integer.MIN_VALUE || out[1] == Integer.MIN_VALUE) {
				out[0] = Integer.MIN_VALUE;
				out[1] = Integer.MIN_VALUE;
			}
			else {
				break;
			}
		}
		return out;
	}
	
	/**
	 * put XA alignments into segments
	 * @param seg_map all segments(output);
	 * @param xa string of "(XA:)...";
	 * @param r removing rRNA;
	 */
	private void putXAstate(HashMap<Integer, ArrayList<Segment>> seg_map, String xa, RemoverRNA r) {
		String[] reads = xa.split(";");
		for (int i = 0; i < reads.length; ++i) {
			boolean p = true;
			String[] cols = reads[i].split(",");
			int start_pos = Integer.parseInt(cols[1]);
			if (start_pos < 0) {
				start_pos = -start_pos;
				p = false;
			}
			int chr_num = ExonInfo.chrSymbolToNum(cols[0]);
			if (chr_num >= 0) {
				int[] cigar = getSegment(cols[2], p);
				int seg_key = cigar[2] * this.read_lenth + cigar[0];
				ArrayList<Segment> seg_list = null;
				if (seg_map.containsKey(seg_key)) {
					seg_list = seg_map.get(seg_key);
				}
				else {
					seg_list = new ArrayList<>();
				}
				Segment seg = new Segment();
				seg.chr_num = chr_num;
				seg.positive = p;
				seg.pos_start = start_pos;
				seg.pos_end = start_pos + cigar[1];
				if (!r.isrRNA(chr_num, seg.pos_start, seg.pos_end)) {
					seg_list.add(seg);
					seg_map.put(seg_key, seg_list);
				}
			}
		}
	}
	
	/**
	 * get all possible clip of this segment, and find back junctions in these clips
	 * @param id id of the these alignments;
	 * @param the_seg key of this segment;
	 * @param dev clipping deviation permitted;
	 * @param segs all segments;
	 * @param seg_list a list to store back junctions;
	 * @param itree a tree to store one site mapping alignments;
	 * @param examin a list for checking pair end;
	 */
	private void getAllClip(String id, int the_seg, int dev, HashMap<Integer, ArrayList<Segment>> segs, ArrayList<Segment> seg_list, ArrayList<IntervalTree<ReadInfo>> itree, ArrayList<Segment> examin) {
		ArrayList<Segment> seg_temp = new ArrayList<>();
		int start = the_seg / this.read_lenth;
		int end = the_seg % this.read_lenth;
		boolean front_seg = false;
		boolean behind_seg = false;
		for (Entry<Integer, ArrayList<Segment>> entry : segs.entrySet()) {
			int key = entry.getKey();
			if (key >= (end-dev) * this.read_lenth && key < (end+dev+1) * this.read_lenth) {
				boolean front_clip = false;
				if (segs.containsKey(the_seg)) {
					ArrayList<Segment> up_seg = segs.get(the_seg);
					ArrayList<Segment> down_seg = entry.getValue();
					for (int i = 0; i < up_seg.size(); ++i) {
						for (int j = 0; j < down_seg.size(); ++j) {
							Segment down = down_seg.get(j);
							Segment up = up_seg.get(i);
							if (up == down) {
								continue;
							}
							if (down.chr_num == up.chr_num && !(up.positive ^ down.positive)) {
								behind_seg = true;
								if (up.positive && down.pos_start <= up.pos_end) {
									Segment seg = new Segment();
									seg.chr_num = up.chr_num;
									seg.pos_start = down.pos_start;
									seg.pos_end = up.pos_end;
									seg_temp.add(seg);
								}
								else if ((!up.positive) && up.pos_start <= down.pos_end) {
									Segment seg = new Segment();
									seg.chr_num = up.chr_num;
									seg.pos_start = up.pos_start;
									seg.pos_end = down.pos_end;
									seg_temp.add(seg);
								}
								else {
									front_clip = true;
								}
							}
						}
					}
					if (!front_clip) {
						seg_list.addAll(seg_temp);
					}
					seg_temp.clear();
				}
			}
			if (start <= key % this.read_lenth + dev && start >= key % this.read_lenth - dev) {
				ArrayList<Segment> up_seg = segs.get(the_seg);
				ArrayList<Segment> down_seg = entry.getValue();
				for (int i = 0; i < up_seg.size(); ++i) {
					for (int j = 0; j < down_seg.size(); ++j) {
						Segment down = down_seg.get(j);
						Segment up = up_seg.get(i);
						if (down.chr_num == up.chr_num && !(up.positive ^ down.positive)) {
							front_seg = true;
						}
					}
				}
			}
		}
		if (!behind_seg && end < this.read_lenth - dev) {
			for (int j = 0; j < segs.get(the_seg).size(); ++j) {
				Segment seg = segs.get(the_seg).get(j);
				boolean front = examin.size() == 0;
				for (int i = 0; i < examin.size(); ++i) {
					if (seg.chr_num == examin.get(i).chr_num && seg.pos_end >= examin.get(i).pos_end) {
						front |= true;
						break;
					}
				}
				if (front) {
					ReadInfo old = itree.get(seg.chr_num).put(seg.pos_start, seg.pos_end, null);
					ReadInfo value = old;
					value.incBehind();
					itree.get(seg.chr_num).put(seg.pos_start, seg.pos_end, value);
				}
			}
		}
		if (!front_seg && the_seg > dev * this.read_lenth) {
			for (int j = 0; j < segs.get(the_seg).size(); ++j) {
				Segment seg = segs.get(the_seg).get(j);
				boolean behind = examin.size() == 0;
				for (int i = 0; i < examin.size(); ++i) {
					if (seg.chr_num == examin.get(i).chr_num && seg.pos_start <= examin.get(i).pos_start) {
						behind |= true;
						break;
					}
				}
				if (behind) {
					ReadInfo old = itree.get(seg.chr_num).put(seg.pos_start, seg.pos_end, null);
					ReadInfo value = old;
					value.incFront();
					itree.get(seg.chr_num).put(seg.pos_start, seg.pos_end, value);
				}
			}
		}
	}
	
	public int getRead_lenth() {
		return read_lenth;
	}

	public void setRead_lenth(int read_lenth) {
		this.read_lenth = read_lenth;
	}
	
}
