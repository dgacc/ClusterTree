package random;

import java.util.Random;

public class MyRandom {
	public static Random r = new Random();
	public static  void  setSeed( int s){
		r.setSeed(s);
	}

}
