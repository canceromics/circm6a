package mapping;

import java.util.HashMap;
import java.util.Iterator;

import genome.IntRegion;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import main.InParam;
import main.Method;

public class Alignment extends IntRegion{

	private String id = null;
	private String seq = null;
	private HashMap<String, IntervalTree<Segment>> segs = null;
	private Alignment pair_end = null;
	private int valid_seg = 0;
	
	private static final int SEG_FORWARD = 1024;
	private static final int MAP_DISTANCE = 50000;
	
	public Alignment(int start, int end, String id) {
		super(start, end);
		this.id = id;
		this.segs = new HashMap<>();
	}
	
	public boolean addSegs(SAMRecord record, int map_q, boolean uniq_map) {
		if (seq == null || record.getReadString().length() > seq.length()) {
			seq = record.getReadString();
		}
		if (!uniq_map || record.getAttribute("XA") == null){
			String key = record.getReferenceName();
//			if (seq == null || record.getReadString().length() > seq.length()) {
//				seq = record.getReadString();
//			}
//			else if (record.getReadString().length() == seq.length() && !seq.equals(record.getReadString())) {
//				System.out.println("debug");
//			}
			putSegment(key, record.getCigarString(), record.getStart(), record.getReadNegativeStrandFlag(), record.getMappingQuality());
			if (InParam.getParams().isXA_mode()) {
				putXAState((String) record.getAttribute("XA"), record.getMappingQuality());
			}
			++valid_seg;
			return true;
		}
		return false;
	}
	
	private void putSegment(String key, String cigar, int start, boolean negative, int quali) {
		int[] cigar_value = cigarToRegion(cigar);
		if (Method.rRNA() != null && Method.rRNA().isrRNA(key, start, start + cigar_value[2] - 1)) {
			MappingMethod.putMappingStat(null, 0, 0, Method.isIP(), quali);
			return;
		}
		Segment s = new Segment(start, start + cigar_value[2] - 1, negative, quali);
		expandRead(cigar_value[0], cigar_value[1]);
		if (!segs.containsKey(key)) {
			segs.put(key, new IntervalTree<>());
		}
		s.addAnotherSeg(segs.get(key).put(cigar_value[0] - SEG_FORWARD, cigar_value[1] - SEG_FORWARD, s));
	}
	
	public static int[] cigarToRegion(String cigar) {
		boolean start_flag = false;
		int start = 1;
		int len = 0;
		int end = 0;
		int seg = 0;
		for (int i = 0; i < cigar.length(); ++i) {
			char c = cigar.charAt(i);
			switch (c) {
			case '0':
			case '1':
			case '2':
			case '3':
			case '4':
			case '5':
			case '6':
			case '7':
			case '8':
			case '9':
				seg = seg * 10 + c - '0';
				break;
				
			case 'M':
				len += seg;
				
			case 'I':
				start_flag = true;
				end += seg;
				seg = 0;
				break;
				
			case 'N':
			case 'D':
				len += seg;
				seg = 0;
				break;
				
			case 'S':
			case 'H':
				start += start_flag ? 0 : seg;
				end += start_flag ? 0 : seg;
				seg = 0;
				break;
				
			default:
				System.err.printf("Warning: %s has undefined char: %c\n", cigar, c);
				seg = 0;
				break;
			}
		}
		return new int[] {start, end, len};
	}
	
	private void expandRead(int start, int end) {
		start = getStart() == 0 ? start : getStart() < start ? getStart() : start;
		end = getEnd() < end ? end : getEnd();
		resetStartAndEnd(start, end);
	}
	
	private void putXAState(String XA, int quali) {
		if (XA == null) {
			return;
		}
		String[] mappings = XA.split(";");
		for (int i = 0; i < mappings.length; ++i) {
			putXAValue(mappings[i], quali);
		}
	}
	
	private void putXAValue(String XA, int quali) {
		String[] cols = XA.split(",");
		if (cols.length >= 3) {
			int start = Integer.parseInt(cols[1]);
			boolean negative = start < 0;
			String key = cols[0];
			putSegment(key, cols[2], negative ? -start : start, negative, quali);
		}
	}
	
	public void handleReads(HashMap<String, IntervalTree<Junction>> junc_table, HashMap<String, IntervalTree<MappingStat>> mapping_tree,
			HashMap<String, IntervalTree<MissingAlignment>> missing_table, boolean search) {
		if (search && valid_seg >= 2) {
			callBackJunction(junc_table, mapping_tree);
		}
		putMappingStat(mapping_tree);
		if (search && seq != null && getLength() + Method.getDev() < seq.length()) {
			callMissingAlignment(missing_table);
		}
		if (pair_end != null && pair_end.valid_seg >= 0) {
			valid_seg = -1;
			pair_end.handleReads(junc_table, mapping_tree, missing_table, search);
		}
		clear();
	}

	private void callBackJunction(HashMap<String, IntervalTree<Junction>> junc_table, HashMap<String, IntervalTree<MappingStat>> mapping_tree) {
		segs.forEach((chr, seg_tree) -> {
			Iterator<Node<Segment>> nodes = seg_tree.overlappers(-SEG_FORWARD, -1);
			while (nodes.hasNext()) {
				Node<Segment> node = nodes.next();
				Iterator<Node<Segment>> next_behinds = seg_tree.overlappers(node.getEnd() + 1 - Method.getDev(), node.getEnd() + 1 + Method.getDev());
				while (next_behinds.hasNext()) {
					Node<Segment> next_behind = next_behinds.next();
					if (next_behind.getStart() >= node.getEnd() + 1 - Method.getDev()) {
						if (!junc_table.containsKey(chr)) {
							if (InParam.getParams().getCirc_bed() != null) {
								continue;
							}
							junc_table.put(chr, new IntervalTree<>());
						}
						if (!mapping_tree.containsKey(chr)) {
							mapping_tree.put(chr, new IntervalTree<>());
						}
						putBackJunction(junc_table.get(chr), mapping_tree.get(chr), callBackJunction(chr, node, next_behind));
					}
				}
			}
		});
	}
	
	private void putBackJunction(IntervalTree<Junction> junc_tree, IntervalTree<MappingStat> mapping_tree, HashMap<Junction, Segment[]> juncs) {
		if (juncs != null) {
			juncs.forEach((junction, segments) -> {
				if (InParam.getParams().getCirc_bed() != null) {
					Node<Junction> junc_node = junc_tree.find(junction.getStart(), junction.getEnd());
					junc_node.getValue().merge(junction);
				}
				else {
					junction.merge(junc_tree.put(junction.getStart(), junction.getEnd(), junction));
				}
				if (!segments[0].isOnceCirc()) {
					MappingMethod.putCircMappingStat(mapping_tree, segments[0].getStart(), segments[0].getEnd(), Method.isIP());
				}
				if (!segments[1].isOnceCirc()) {
					MappingMethod.putCircMappingStat(mapping_tree, segments[1].getStart(), segments[1].getEnd(), Method.isIP());
				}
			});
		}
	}

	private HashMap<Junction, Segment[]> callBackJunction(String key, Node<Segment> front, Node<Segment> behind){
		HashMap<Junction, Segment[]> out = new HashMap<>();
		Segment seg = front.getValue();
		while (seg != null) {
			Segment behind_seg = behind.getValue();
			while (behind_seg != null) {
				if (behind_seg != seg && !seg.isNegative() ^ behind_seg.isNegative() 
						&& (seg.isEnoughQuality() || behind_seg.isEnoughQuality())) {
					if (seg.getEnd() > behind_seg.getStart()) {
						if (!hasPairFront(key, behind_seg.getStart()) && !hasPairBehind(key, seg.getEnd())
								&& Method.withIn(seg.getEnd() - behind_seg.getStart() + 1, InParam.getParams().getCirc_length(), InParam.getParams().getCirc_max_length())) {
							Junction junction = getBackJunction(key, front, behind, seg, behind_seg);
							if (junction != null) {
								out.put(junction, new Segment[] {seg, behind_seg});
							}
						}
					}
					else if (seg.getEnd() + MAP_DISTANCE >= behind_seg.getStart()){
						return null;
					}
				}
				behind_seg = behind_seg.getAnotherSeg();
			}
			seg = seg.getAnotherSeg();
		}
		return out;
	}
	
	private boolean hasPairFront(String key, int start) {
		return pair_end != null && pair_end.segs.containsKey(key) && pair_end.segs.get(key).overlappers(Math.max(1, start - MAP_DISTANCE), start - Method.getDev() - 1).hasNext();
	}

	private boolean hasPairBehind(String key, int end) {
		return pair_end != null && pair_end.segs.containsKey(key) && pair_end.segs.get(key).overlappers(end + Method.getDev() + 1, end + MAP_DISTANCE).hasNext();
	}
	
	private Junction getBackJunction(String chr, Node<?> first, Node<?> last, Segment head, Segment tail) {
		int[] adjust = null;
		int adjust_start = Math.min(first.getEnd(), last.getStart() - 1);
		int adjust_end = Math.max(first.getEnd(), last.getStart() - 1);
		for (int i = adjust_start; i <= adjust_end; i++) {
				adjust = Method.getInstance().adjustToClip(chr, tail.getStart() + i + 1 - last.getStart(), head.getEnd() + i - first.getEnd(), 
						seq.substring(i + SEG_FORWARD, last.getEnd() + SEG_FORWARD), seq.substring(first.getStart() - 1 + SEG_FORWARD, i + SEG_FORWARD));
			if (adjust != null) {
//				if (adjust[0] < 1 || adjust[1] < 0 || adjust[0] > adjust[1]) {
//					System.out.println("Debug:");
//					System.out.printf("%s\nFirst:\n%d\t%d\t%d\t%d\nLast:\n%d\t%d\t%d\t%d\n",
//							this.id, first.getStart(), first.getEnd(), head.getStart(), head.getEnd(),
//							last.getStart(), last.getEnd(), tail.getStart(), tail.getEnd());
//				}
				if ((head.isEnoughQuality() || Method.getInstance().isCircClip(Method.getChrSeq(chr, adjust[0] - 1, adjust[1]), 
						Method.getChrSeq(chr, Math.max(0, adjust[0] - Math.max(50000, adjust[1] - adjust[0] + Method.getDev())), adjust[0] - 1), 
						seq.substring(first.getStart() - 1 + SEG_FORWARD, i + SEG_FORWARD)))
						&& (tail.isEnoughQuality() || Method.getInstance().isCircClip(Method.getChrSeq(chr, adjust[0] - 1, adjust[1]),
								Method.getChrSeq(chr, adjust[1], Math.min(Method.getChrLength(chr), adjust[1] + Math.max(50000, adjust[1] - adjust[0] + Method.getDev()))), 
								seq.substring(i + SEG_FORWARD, last.getEnd() + SEG_FORWARD)))) {
					break;
				}
				else {
					return null;
				}
			}
		}
		if (adjust != null && Method.withIn(adjust[1] - adjust[0] + 1, InParam.getParams().getCirc_length(), InParam.getParams().getCirc_max_length())) {
			Junction junc = new Junction(adjust[0], adjust[1]);
			junc.addID(id, Method.isIP());
			return junc;
		}
		return null;
	}
	
	private void putMappingStat(HashMap<String, IntervalTree<MappingStat>> mapping_tree) {
		segs.forEach((chr, seg_tree) -> {
			Iterator<Node<Segment>> nodes = seg_tree.overlappers(-SEG_FORWARD, -1);
			while (nodes.hasNext()) {
				Node<Segment> node = nodes.next();
				if (!mapping_tree.containsKey(chr)) {
					mapping_tree.put(chr, new IntervalTree<>());
				}
				for (Segment s = node.getValue();s != null; s = s.getAnotherSeg()) {
					MappingMethod.putMappingStat(mapping_tree.get(chr), s.getStart(), s.getEnd(), Method.isIP(), s.getQuality());
				}
			}
		});
	}

	private void callMissingAlignment(HashMap<String, IntervalTree<MissingAlignment>> missing_table) {
		if (getStart() >= Method.getDev()) {
			segs.forEach((chr, seg_tree) -> {
				if (!missing_table.containsKey(chr)) {
					missing_table.put(chr, new IntervalTree<>());
				}
				Iterator<Node<Segment>> nodes = seg_tree.overlappers(getStart() - SEG_FORWARD, getStart() - SEG_FORWARD + Method.getDev());
				while (nodes.hasNext()) {
					Node<Segment> node = nodes.next();
					for (Segment s = node.getValue(); s != null; s = s.getAnotherSeg()) {
						MissingAlignment mis = new MissingAlignment(true, Method.isIP(), seq.substring(0, node.getStart() + SEG_FORWARD - 1), s);
						mis.addAnotherMissingAlignment(missing_table.get(chr).put(s.getStart() - Method.getDev(), s.getStart() + Method.getDev(), mis));
					}
				}
			});
			
		}
		if (getEnd() <= seq.length() - Method.getDev()) {
			segs.forEach((chr, seg_tree) -> {
				if (!missing_table.containsKey(chr)) {
					missing_table.put(chr, new IntervalTree<>());
				}
				Iterator<Node<Segment>> nodes = seg_tree.overlappers(getEnd() - SEG_FORWARD - Method.getDev(), getEnd() - SEG_FORWARD);
				while (nodes.hasNext()) {
					Node<Segment> node = nodes.next();
					for (Segment s = node.getValue(); s != null; s = s.getAnotherSeg()) {
						MissingAlignment mis = new MissingAlignment(false, Method.isIP(), seq.substring(node.getEnd() + SEG_FORWARD, seq.length()), s);
						mis.addAnotherMissingAlignment(missing_table.get(chr).put(s.getEnd() - Method.getDev(), s.getEnd() + Method.getDev(), mis));
					}
				}
			});
		}
	}
	
	public String getID() {
		return id;
	}
	
	public Alignment getPairRead() {
		return pair_end;
	}
	
	public void setPairRead(Alignment pair_end) {
		this.pair_end = pair_end;
	}
	
	public void clear() {
		pair_end = null;
		segs.clear();
		segs = null;
		id = null;
		seq = null;
	}
}
