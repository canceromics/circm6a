package genome;

public class Exon extends IntRegion{

	private Transcript script = null;
	
	public Exon(Transcript script, int start, int end) {
		super(start, end);
		this.script = script;
	}

	public String getChr() {
		return script==null? null : script.getChr();
	}
	
	public Gene getGene() {
		return script==null? null : script.getGene();
	}
	
	public Transcript getScript() {
		return script;
	}
	
	public int fixToExonStart(int target, int dev) {
		return Math.abs(getStart() - target) <= dev ? getStart() : -128;
	}
	
	public int fixToExonEnd(int target, int dev) {
		return Math.abs(getEnd() - target) <= dev ? getEnd() : -128;
	}
}
