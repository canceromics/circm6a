package main;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;

import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.IntervalTree.Node;

public class RemoverRNA {
	
	private HashMap<String, IntervalTree<Integer>> rRNA_trees =null;
	
	/**
	 * using the bed file to build filting tree, no tree if null
	 * @param bed rRNA bed file path
	 */
	public RemoverRNA(String bed) {
		if (bed != null) {
			this.rRNA_trees = new HashMap<>();
			this.creatTree(bed);
		}
	}
	
	/**
	 * determine whether this region is a rRNA region
	 * @param chr region chromosome
	 * @param start region start
	 * @param end region end
	 * @return true if no less than a half drop in rRNA region
	 */
	public boolean isrRNA(String chr, int start, int end) {
		boolean out = false;
		if (this.rRNA_trees != null) {
			IntervalTree<Integer> tree = this.rRNA_trees.get(chr);
			if (tree == null){
				return out;
			}
			Iterator<Node<Integer>> nodes = tree.overlappers(start, end);
			while (nodes.hasNext()) {
				Node<Integer> node = nodes.next();
				if ((Math.min(end, node.getEnd()) - Math.max(start, node.getStart())) * 2 >= end - start) {
					out = true;
					break;
				}
			}
		}
		return out;
	}
	
	/**
	 * creat filting tree from bed file
	 * @param file file path
	 */
	public void creatTree(String file) {
		InputStreamReader isr = null;
		BufferedReader reader = null;
		try {
			if (file.equals("/rRNA.bed")) {
				InputStream is = Object.class.getResourceAsStream(file);
				isr = new InputStreamReader(is, "UTF-8");
			}
			else {
				isr = new InputStreamReader(new FileInputStream(file), "UTF-8");
			}
			reader = new BufferedReader(isr);
			String line = null;
			while ((line = reader.readLine()) != null) {
				String[] cols = line.split("\t");
				String chr = cols[0];
				if (this.rRNA_trees.containsKey(chr)){
					IntervalTree<Integer> temp_T = new IntervalTree<>();
					temp_T.setSentinel(null);
					this.rRNA_trees.put(chr, temp_T);
				}
				IntervalTree<Integer> tree = this.rRNA_trees.get(chr);
				int start = Integer.parseInt(cols[1]);
				int end = Integer.parseInt(cols[2]);
				tree.put(start, end, null);
			}
		}
		catch(IOException e){
			e.printStackTrace();
		}
		finally{
			if (reader != null){
				try{
					reader.close();
				}
				catch(IOException e1){
				}
			}
		}
	}
}
