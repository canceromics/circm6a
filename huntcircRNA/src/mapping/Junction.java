package mapping;

import java.util.HashSet;

import genome.Gene;
import genome.IntRegion;
import main.InParam;

public class Junction extends IntRegion{

	private HashSet<String> ids = new HashSet<>(4);
	private HashSet<String> ids_ip = new HashSet<>(4);
	private int juncs_ip = 0;
	private int juncs = 0;
	private int loose_ip = 0;
	private int loose = 0;
	private Gene gene = null;
	private boolean fix_start = false;
	private boolean fix_end = false;
	private boolean visited = false;
	
	public Junction(int start, int end) {
		super(start, end);
	}
	
	public boolean addID(String id, boolean ip_flag) {
		return ip_flag ? ids_ip.add(id) : ids.add(id);
	}
	
	public int incJunc(boolean ip_flag) {
		return ip_flag ? ++juncs_ip : ++juncs;
	}
	
	public int incLoose(boolean ip_flag) {
		return ip_flag ? ++loose_ip : ++loose;
	}
	
	public HashSet<String> getIDs(){
		return ids;
	}
	
	public HashSet<String> getIPIDs(){
		return ids_ip;
	}
	
	public String getIDString() {
		StringBuilder sb = new StringBuilder();
		for (String id : ids) {
			sb.append(id);
			sb.append(',');
		}
		if (sb.length() > 0) {
			sb.setLength(sb.length() - 1);
			sb.append(';');
		}
		for (String id : ids_ip) {
			sb.append(id);
			sb.append(',');
		}
		if (sb.length() > 0) {
			sb.setLength(sb.length() - 1);
		}
		return sb.toString();
	}
	
	public void fixStart() {
		fix_start = true;
	}
	
	public boolean isStartFixed() {
		return fix_start;
	}
	
	public void fixEnd() {
		fix_end = true;
	}
	
	public boolean isEndFixed() {
		return fix_end;
	}
	
	public boolean isAllFixed() {
		return fix_start && fix_end;
	}
	
	public void visit() {
		visited = true;
	}
	
	public boolean isVisited() {
		return visited;
	}
	
	public void setGene(Gene gene) {
		this.gene = gene;
	}
	
	public void resetGene(Gene gene) {
		this.gene = this.gene == null || (gene.getLength() > this.gene.getLength()) ? gene : this.gene;
	}
	
	public Gene getGene() {
		return gene;
	}
	
	public int getJuncNum() {
		return juncs;
	}
	
	public int getJuncNumIP() {
		return juncs_ip;
	}
	
	public int getLooseIP() {
		return loose_ip;
	}
	
	public int getLoose() {
		return loose;
	}
	
	public int getConfidence(int t1, int t2) {
		return (isIPpropor() && ids_ip.size() + juncs_ip >= t1) ? 1 : (ids_ip.size() + juncs_ip > 0 ? 0 : -1); 
	}
	
	private boolean isIPpropor() {
		if (ids_ip.size() + juncs_ip == 0) {
			return false;
		}
		if (ids.size() + juncs == 0) {
			return true;
		}
		return (double) (ids_ip.size() + juncs_ip) / (ids.size() + juncs) * MappingStat.getMappingInputReads()
				/ MappingStat.getMappingIpReads() >= InParam.getParams().getPropor();
	}
	
	public void merge(Junction junction) {
		if (junction != null) {
			this.ids.addAll(junction.ids);
			this.ids_ip.addAll(junction.ids_ip);
			this.juncs += junction.juncs;
			this.juncs_ip += junction.juncs_ip;
		}
	}
	
	public boolean checkBounderAndGene(boolean gene_mode) {
		return fix_start && fix_end && !(gene_mode && gene == null);
	}
	
	public boolean reachOutputStandard(int single) {
		return ids.size() + ids_ip.size() + juncs + juncs_ip + single >= InParam.getParams().getSup_reads();
	}
	
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(juncs);
		sb.append('\t');
		if (InParam.getParams().getIp_file() != null) {
			sb.append(juncs_ip);
			sb.append('\t');
		}
		sb.append(getIDString());
		return sb.toString();
	}
}
