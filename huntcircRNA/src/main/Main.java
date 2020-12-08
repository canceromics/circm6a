package main;

public class Main {
	
	public static void main (String[] args){		
		InParam.getParams().SetParams(args);
		if (InParam.getParams().completeArgs()) {
			Method.printNow("Start at");
			Method.run();
			
			Method.printNow("Finish at");
		}
	}
	
}