package mapping;

import java.util.HashMap;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class MappingMethod {

	public static void putMappingStat(IntervalTree<MappingStat> mapping_tree, int start, int end, boolean ip, int mapq) {
		if (mapping_tree != null) {
			MappingStat ms = mapping_tree.put(start, end, null);
			if (ms == null) {
				ms = new MappingStat();
			}
			ms.incReads(ip, mapq);
			mapping_tree.put(start, end, ms);
		}
		MappingStat.incReadNum(ip, mapq);
	}
	
	
	public static void putCircMappingStat(IntervalTree<MappingStat> mapping_tree, int start, int end, boolean ip) {
		MappingStat ms = mapping_tree.put(start, end, null);
		if (ms == null) {
			ms = new MappingStat();
		}
		ms.incCircReads(ip);
		mapping_tree.put(start, end, ms);
	}
	
	public static void resetMappingStat(HashMap<String, IntervalTree<MappingStat>> mapping_table, boolean ip) {
		MappingStat.resetReadNum(ip);
		if (mapping_table != null) {
			mapping_table.forEach((chr, tree) -> {
				resetMappingTree(tree, ip);
			});
		}
	}
	
	public static void resetMappingTree(IntervalTree<MappingStat> mapping_tree, boolean ip) {
		if (mapping_tree != null) {
			for (Node<MappingStat> ms_node : mapping_tree) {
				resetNode(ms_node.getValue(), ip);
			}
		}
		MappingStat.resetReadNum(ip);
	}
	
	public static void resetNode(MappingStat ms, boolean ip) {
		if (ms != null) {
			ms.resetReads(ip);
		}
	}
}
