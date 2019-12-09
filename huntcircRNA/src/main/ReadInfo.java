package main;

public class ReadInfo {

	private int no_xa = 0;
//	private int xa = 0;
	private int front = 0;
	private int behind = 0;
	private int circ_front = 0;
	private int circ_behind = 0;
	private int linear = 0;
	
	public ReadInfo() {
	}
	
	public int getNo_xa() {
		return no_xa;
	}
	
	public void incNo_xa() {
		++this.linear;
		++this.no_xa;
	}
	
	public void setNo_xa(int no_xa) {
		this.no_xa = no_xa;
	}
	
	public int getFront() {
		return front;
	}

	public void incFront() {
		++this.front;
	}
	
	public void setFront(int front) {
		this.front = front;
	}

	public int getBehind() {
		return behind;
	}

	public void incBehind() {
		++this.behind;
	}
	
	public void setBehind(int behind) {
		this.behind = behind;
	}
	
	public int getCirc_front() {
		return circ_front;
	}

	public int getCirc_behind() {
		return circ_behind;
	}

	public void incCirc_front() {
		++this.circ_front;
	}
	
	public void incCirc_behind() {
		++this.circ_behind;
	}
	
	public int getLinear() {
		return linear;
	}

	public void calLinear(boolean start_flag) {
		if (start_flag) {
			this.linear = Math.min(this.linear, this.no_xa - this.circ_front - this.front);
		}
		else {
			this.linear = Math.min(this.linear, this.no_xa - this.circ_behind - this.behind);
		}
		if (this.linear < 0) {
			System.out.println("Linear Reads Warning Debug");
			this.linear = 0;
		}
	}
	
	public void setLinear(int linear) {
		this.linear = linear;
	}
	
	
}
