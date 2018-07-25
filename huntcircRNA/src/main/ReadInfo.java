package main;

import java.util.ArrayList;

public class ReadInfo {

	private int no_xa = 0;
//	private int xa = 0;
	private int front = 0;
	private int behind = 0;
	private int circ = 0;
	private int linear = 0;
	private ArrayList<String> no_front = null;
	private ArrayList<String> no_behind= null;
	
	public ReadInfo() {
	}
	
	public ReadInfo(boolean list) {
		if (list) {
			this.no_front = new ArrayList<>();
			this.no_behind = new ArrayList<>();
		}
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
	
//	public int getXa() {
//		return xa;
//	}
//	
//	public void incXa() {
//		++this.xa;
//	}
//	
//	public void setXa(int xa) {
//		this.xa = xa;
//	}
	
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

	public int getCirc() {
		return circ;
	}

	public void incCirc() {
		++this.circ;
	}
	
	public int getLinear() {
		return linear;
	}

	public void calLinear(boolean start_flag) {
		if (start_flag) {
			this.linear = Math.min(this.linear, this.no_xa - this.circ - this.front);
		}
		else {
			this.linear = Math.min(this.linear, this.no_xa - this.circ - this.behind);
		}
	}
	
	public void setLinear(int linear) {
		this.linear = linear;
	}

	public ArrayList<String> getNo_front() {
		return no_front;
	}
	
	public void addNo_front(String id) {
		this.no_front.add(id);
	}
	
	public void setNo_front(ArrayList<String> no_front) {
		this.no_front = no_front;
	}
	
	public ArrayList<String> getNo_behind() {
		return no_behind;
	}
	
	public void addNo_behind(String id) {
		this.no_behind.add(id);
	}
	
	public void setNo_behind(ArrayList<String> no_behind) {
		this.no_behind = no_behind;
	}
	
	
}
