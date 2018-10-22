package structures;

import java.util.LinkedList;

public class Cycles {
	private LinkedList<Integer> cycle= new LinkedList<Integer>();
	private int maxElement;
	
	public LinkedList<Integer> getCycle() {
		return cycle;
	}

	public void setCycle(LinkedList<Integer> cycle) {
		this.cycle = cycle;
	}

	public int getMaxElement() {
		return maxElement;
	}

	public void setMaxElement(int maxElement) {
		this.maxElement = maxElement;
	}

	public LinkedList<Integer>  sortElement(int maxElement){
		 int index  = cycle.indexOf(maxElement);
		 int length = cycle.size();
		    for( int j = length - 1; j > index; j--){
		    	if(cycle.get(j) < maxElement){
		    		cycle.addFirst(cycle.get(j));
		    		cycle.remove(length);
		    	}
		    }
		    return  cycle;
	}
	

}
