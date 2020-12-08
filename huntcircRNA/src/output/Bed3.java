package output;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;

import genome.IntRegion;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;
import main.Method;

public class Bed3 extends IntRegion{
	
	public final static String HEADER = "#Chr\tStart\tEnd";
	
	private static ArrayList<String> chr_in_order = null;
	private static char sep = '\t';
	
	private String chr = null;
	
	public Bed3(String chr, int start, int end) {
		super(start, end);
		this.chr = chr;
	}

	public String getChr() {
		return chr;
	}

	public static char getSep() {
		return sep;
	}

	public static void setSep(char sep) {
		Bed3.sep = sep;
	}

	public static ArrayList<String> getChrsInOrder() {
		return chr_in_order;
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(chr);
		sb.append(sep);
		sb.append(getStart());
		sb.append(sep);
		sb.append(getEnd());
		return sb.toString();
	}
	
	public static <T extends Collection<String>> boolean setChrOrder(T chrs) {
		if (chr_in_order == null) {
			chr_in_order = new ArrayList<>();
			chr_in_order.addAll(chrs);
			Collections.sort(chr_in_order, new Comparator<String>() {

				@Override
				public int compare(String o1, String o2) {
					if (o1.length() != o2.length()) {
						int l = Math.min(o1.length(), o2.length()) - 1;
						if (o1.charAt(l) <= '9' && o1.charAt(l) >= '0'
								&& o2.charAt(l) <= '9' && o2.charAt(l) >= '0') {
							return o1.length() - o2.length();
						}
					}
					return o1.compareTo(o2);
				}
				
			});
		}
		return false;
	}
	
	public static <T> void writeBedInOrder(BufferedWriter writer, HashMap<String, IntervalTree<T>> out) throws IOException {
		if (chr_in_order == null) {
			Method.writeTreeMap(writer, out);
		}
		else {
			for (String key : chr_in_order) {
				if (out.get(key) == null) {
					continue;
				}
				Iterator<Node<T>> nodes = out.get(key).overlappers(Integer.MIN_VALUE, Integer.MAX_VALUE);
				while (nodes.hasNext()) {
					Node<T> node = nodes.next();
					writer.write(key);
					writer.write(Bed3.getSep());
					writer.write(node.getStart());
					writer.write(Bed3.getSep());
					writer.write(node.getEnd());
					writer.write(Bed3.getSep());
					writer.write(node.getValue().toString());
					writer.newLine();
				}
				writer.flush();
			}
		}
	}
}
